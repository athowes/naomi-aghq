tmb_input <- tmb_inputs
DLL <- "naomi_simple"
inner_verbose <- FALSE
progress <- NULL
map <- NULL

warning("For the aghq, k must be set to 1 as thing stand currently.")
stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)
quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = 1)

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
theta <- as.numeric(quad$modesandhessians[theta_names])

.f <- function(i) {
  tmb_input_i$data$i <- i
  obj_i <- local_make_tmb_obj(tmb_input_i$data, tmb_input_i$par, calc_outputs = 0L, inner_verbose, progress, DLL = "naomi_simple_x_index")
  random_i <- obj_i$env$random
  mode_i <- quad$modesandhessians[["mode"]][[1]][-i]
  obj_i$env$last.par[random_i] <- mode_i

  out <- data.frame(
    index = i,
    par = x_names[i],
    x = spline_nodes(quad$modesandhessians, i = i, k = 5)
  )

  .g <- function(x) {
    lp <- as.numeric(- obj_i$fn(c(x, theta))) + log(quad$normalized_posterior$nodesandweights$weights)
    return(lp - quad$normalized_posterior$lognormconst)
  }

  out$lp <- purrr::map_dbl(out$x, .g)
  return(out)
}

out <- purrr::map(.x = 1:10, .f = .f)
out <- dplyr::bind_rows(out)
return(out)
