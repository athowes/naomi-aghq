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


