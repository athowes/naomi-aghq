#' Develop a version of adam Laplace marginals method where we use k = 2 but we
#' fix some of the hyperparameters to their EB values using the map option

#' Start of if(adam) loop

#' First find the fitted values of the hyperparameters then use map to fix them
tmb <- out
tmb$fit$par.fixed

#' Suppose that we fix everything but:
#' * `logit_phi_rho_x`
#' * `logit_phi_rho_xs`
#'
#' These both have relatively large uncertainty

#' This would give us 3^2 = 9 points
n_hyper <- 2
k <- 3
k^n_hyper

theta_names <- unique(names(tmb$fit$par))
fixed_theta_names <- setdiff(theta_names, c("logit_phi_rho_x", "logit_phi_rho_xs"))
map_fixed_theta <- list()
param_fixed_theta <- tmb_inputs_simple$par_init

for(theta in fixed_theta_names) {
  map_fixed_theta[[theta]] <- rep(factor(NA), length(tmb_inputs_simple$par_init[[theta]]))
  param_fixed_theta[[theta]] <- tmb$fit$par.fixed[[theta]]
}

#' Start expose fit_aghq(tmb_inputs_simple, k = k, basegrid = sparse_grid, control = control)
tmb_input <- tmb_inputs_simple
k <- k
basegrid <- NULL
control <- default_control_tmb()
inner_verbose <- FALSE
progress <- NULL
map <- map_fixed_theta
DLL <- "naomi_simple"

tmb_input$par_init <- param_fixed_theta

obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)

#' Start expose aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, ...)
ff <- obj
startingvalue <- obj$par
transformation <- aghq::default_transformation()
optresults <- NULL
basegrid <- basegrid
control <- control

validate_control(control, type = "tmb")
validate_transformation(transformation)

transformation <- make_transformation(transformation)

thetanames <- NULL

if (exists("par", ff)) thetanames <- make.unique(names(ff$par), sep = "")

if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr, theta, method = "Richardson")
}

quad <- aghq(
  ff = ff, k = k, transformation = transformation,
  startingvalue = startingvalue, optresults = optresults,
  basegrid = basegrid, control = control
)

distinctthetas <- quad$normalized_posterior$nodesandweights[, grep("theta", colnames(quad$normalized_posterior$nodesandweights))]

if (!is.data.frame(distinctthetas)) distinctthetas <- data.frame(theta1 = distinctthetas)
modesandhessians <- distinctthetas

if (is.null(thetanames)) {
  thetanames <- colnames(distinctthetas)
} else {
  colnames(modesandhessians) <- thetanames
  colnames(quad$normalized_posterior$nodesandweights)[grep("theta", colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
}

modesandhessians$mode <- vector(mode = "list", length = nrow(distinctthetas))
modesandhessians$H <- vector(mode = "list", length = nrow(distinctthetas))

for (i in 1:nrow(distinctthetas)) {
  theta <- as.numeric(modesandhessians[i, thetanames])
  ff$fn(theta)
  mm <- ff$env$last.par
  modesandhessians[i, "mode"] <- list(list(mm[ff$env$random]))
  H <- ff$env$spHess(mm, random = TRUE)
  H <- rlang::duplicate(H)
  modesandhessians[i, "H"] <- list(list(H))
}

modesandhessians$weights <- quad$normalized_posterior$nodesandweights$weights

quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")
quad$obj <- obj

#' Now we've fit the AGHQ we can get to the Laplace marginals (again, over and over again)
modesandhessians <- quad$modesandhessians

random <- quad$obj$env$random
x_names <- names(quad$obj$env$par[random])

#' Replace the $par_init of tmb_inputs_simple with a concatenated version
#' which has all the latent field parameters compressed into one x
tmb_inputs_simple_i <- tmb_inputs_simple
tmb_inputs_simple_i$par_init <- param_fixed_theta
x_lengths <- lengths(tmb_inputs_simple_i$par_init[unique(x_names)])
tmb_inputs_simple_i$par_init$x_minus_i <- rep(0, sum(x_lengths) - 1)
tmb_inputs_simple_i$par_init$x_i <- 0
tmb_inputs_simple_i$par_init[unique(x_names)] <- NULL

#' Pass in the following as data
#' * `x_lengths`: length of each subvector within x
#' * `x_starts`: start index (remember zero indexing) of each subvector within x
#' * `i`:
x_starts <- cumsum(x_lengths) - x_lengths
tmb_inputs_simple_i$data$x_lengths <- as.numeric(x_lengths)
tmb_inputs_simple_i$data$x_starts <- as.numeric(x_starts)
tmb_inputs_simple_i$data$i <- 7 #' "beta_anc_rho"

#' Change the DLL to the x_index version
DLL <- "naomi_simple_x_index"
compile(paste0(DLL, ".cpp"))
dyn.load(dynlib(DLL))

data <- tmb_inputs_simple_i$data
par <- tmb_inputs_simple_i$par_init
calc_outputs <- 0L
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL
map <- map_fixed_theta

data$calc_outputs <- as.integer(calc_outputs)

obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = DLL,
  silent = !inner_verbose,
  random = "x_minus_i",
  map = map,
  random.start = expression(last.par[random])
)

if (!is.null(progress)) {
  obj$fn <- naomi:::report_progress(obj$fn, progress)
}

#' The set of input values that we'd like to calculate the log-probability at
theta_mode_location <- which.max(quad$normalized_posterior$nodesandweights$logpost_normalized)
modeandhessian <- modesandhessians[theta_mode_location, ]
gg <- create_approx_grid(modeandhessian, tmb_inputs_simple_i$data$i, k = 7)
nodes <- mvQuad::getNodes(gg)

#' Check the order of parameters in obj: there should just be three `x_i`, `logit_phi_rho_x` and `logit_phi_rho_xs`
obj$par

laplace_marginal <- function(x) {
  lp <- vector(mode = "numeric", length = nrow(modesandhessians))
  random <- obj$env$random   #' Indices of random effects

  for(z in 1:nrow(modesandhessians)) {
    theta <- as.numeric(modesandhessians[z, thetanames])
    mode <- modesandhessians[z, "mode"][[1]][-tmb_inputs_simple_i$data$i]
    obj$env$last.par[random] <- mode
    lp[z] <- as.numeric(- obj$fn(c(x, theta)))
  }

  logSumExpWeights(lp, w = modesandhessians$weights) - quad$normalized_posterior$lognormconst
}

lps <- vector(mode = "numeric", length = length(nodes))
starts <- vector(mode = "numeric", length = length(nodes))
ends <- vector(mode = "numeric", length = length(nodes))
times <- vector(mode = "numeric", length = length(nodes))

for(i in seq_along(nodes)) {
  starts[i] <- Sys.time()
  lps[i] <- laplace_marginal(nodes[i])
  ends[i] <- Sys.time()
  times[i] <- ends[i] - starts[i]
}

pdf("beta_anc_rho-map-marginal.pdf", h = 4, w = 6.25)

plot_marginal_spline(nodes, lps)

dev.off()

#' We could also use the normalising constant de jour calculated fresh
#' I am adding back on quad$normalized_posterior$lognormconst
lognormconst <- logSumExpWeights(lps + quad$normalized_posterior$lognormconst, mvQuad::getWeights(gg))
abs(lognormconst - quad$normalized_posterior$lognormconst)

#' Looks to be getting the right answer to me!
