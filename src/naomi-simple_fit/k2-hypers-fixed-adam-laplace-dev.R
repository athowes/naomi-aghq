#' Develop a version of adam Laplace marginals method where we use k = 2 but we
#' fix some of the hyperparameters to their EB values using the map option or similar

#' Start of if(adam) loop

#' First find the fitted values of the hyperparameters then use map to fix them
tmb <- out
tmb$fit$par.fixed

#' Suppose that we fix everything but:
#' * `logit_phi_rho_x`
#' * `logit_phi_rho_xs`
#'
#' These both have relatively large uncertainty -- see adaptive-k-dev.R

#' This would give us 5^2 = 25 points
n_hyper <- 2
k <- 5

theta_names <- unique(names(tmb$fit$par))
fixed_theta_names <- setdiff(theta_names, c("logit_phi_rho_x", "logit_phi_rho_xs"))
map_fixed_theta <- list()
param_fixed_theta <- tmb_inputs$par_init

for(theta in fixed_theta_names) {
  map_fixed_theta[[theta]] <- rep(factor(NA), length(tmb_inputs$par_init[[theta]]))
  param_fixed_theta[[theta]] <- tmb$fit$par.fixed[[theta]]
}

#' Start expose fit_aghq(tmb_inputs, k = k, basegrid = sparse_grid, control = control)
tmb_input <- tmb_inputs
k <- 5
basegrid <- NULL
control <- default_control_tmb()
inner_verbose <- FALSE
progress <- NULL
map <- map_fixed_theta
DLL <- "naomi_simple"

stopifnot(inherits(tmb_input, "naomi_tmb_input"))

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

#' Replace the $par_init of tmb_inputs_simple with a concatendated version
#' which has all the latent field parameters compressed into one x
tmb_inputs_simple_i <- tmb_inputs_simple
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
map <- NULL

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

#' The set of input values that we'd like to calculate the log-probability at...
