#' (This is the TMB fit)
hypers <- names(fit$par)
length(hypers)

sd_out <- TMB::sdreport(obj = fit$obj, par.fixed = fit$par, getJointPrecision = TRUE)

base_grid <- sd_levels_ghe_grid(
  dim = length(hypers),
  level = c(1, 2),
  cut_off = c(0, 1.9),
  sd = sqrt(diag(sd_out$cov.fixed))
)

prod(base_grid$level)

#' Expose fit_adam_grid(tmb_inputs, basegrid = base_grid)
tmb_input <- tmb_inputs
basegrid <- base_grid
inner_verbose <- FALSE
progress <- NULL
map <- NULL
DLL <- "naomi_simple"

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)
# Can optresults be passed in here from previous TMB fit?
quad <- aghq::marginal_laplace_tmb(obj, basegrid = basegrid, startingvalue = obj$par)
quad$obj <- obj

random <- obj$env$random
x_names <- names(obj$env$par[random])
x_lengths <- lengths(tmb_input$par_init[unique(x_names)])
x_starts <- cumsum(x_lengths) - x_lengths

tmb_input_i <- tmb_input
tmb_input_i$par_init$x_minus_i <- rep(0, sum(x_lengths) - 1)
tmb_input_i$par_init$x_i <- 0
tmb_input_i$par_init[unique(x_names)] <- NULL
tmb_input_i$data$x_lengths <- x_lengths
tmb_input_i$data$x_starts <- x_starts

theta_names <- make.unique(names(obj$par), sep = "")

.f <- function(i) {
  tmb_input_i$data$i <- i
  obj_i <- local_make_tmb_obj(tmb_input_i$data, tmb_input_i$par, calc_outputs = 0L, inner_verbose, progress, DLL = "naomi_simple_x_index")
  random_i <- obj_i$env$random
  mode_i <- quad$modesandhessians[["mode"]][[1]][-i]
  gg <- create_approx_grid(quad$modesandhessians, i = i, k = 5)
  out <- data.frame(index = i, par = x_names[i], x = mvQuad::getNodes(gg), w = mvQuad::getWeights(gg))

  .g <- function(x) {
    lp <- vector(mode = "numeric", length = nrow(quad$modesandhessians))

    for(z in 1:nrow(quad$modesandhessians)) {
      theta <- as.numeric(quad$modesandhessians[z, theta_names])
      obj_i$env$last.par[random_i] <- quad$modesandhessians[z, "mode"][[1]][-tmb_input_i$data$i]
      lp[z] <- as.numeric(- obj_i$fn(c(x, theta)))
    }

    return(logSumExpWeights(lp, w = quad$normalized_posterior$nodesandweights$weights))
  }

  out$lp <- purrr::map_dbl(out$x, .g) # Calculate log-posterior
  lognormconst <- logSumExpWeights(out$lp, out$w) # Calculate log normalising constant
  out$lp_normalised <- out$lp - lognormconst # Calculate normalised log-posterior

  return(out)
}

d <- length(random)
d <- 10 # Test how long it takes for 10

start <- Sys.time()
laplace_marginals <- purrr::map(.x = 1:d, .f = .f)
laplace_marginals <- dplyr::bind_rows(laplace_marginals)
end <- Sys.time()

end - start

out <- list("quad" = quad, "laplace_marginals" = laplace_marginals)

laplace_marginals
