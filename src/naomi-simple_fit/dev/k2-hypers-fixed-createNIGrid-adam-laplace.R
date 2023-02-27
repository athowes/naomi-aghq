#' Develop a version of adam Laplace marginals method where we use k = 2 but we
#' fix some of the hyperparameters to their EB values using a k varied NIGrid

#' Start of if(adam) loop

#' (This is the TMB fit)
hypers <- names(fit$par)
length(hypers)

sd_out <- sdreport(obj = fit$obj, par.fixed = fit$par, getJointPrecision = TRUE)

sd_levels_ghe_grid <- function(dim, level, cut_off, sd) {
  stopifnot(length(level) == length(cut_off))
  stopifnot(dim == length(sd))
  levels <- vector(mode = "numeric", length = dim)
  for(i in seq_along(cut_off)) levels[sd > cut_off[i]] <- level[i]
  grid <- mvQuad::createNIGrid(dim = dim, "GHe", level = levels)
  grid
}

base_grid <- sd_levels_ghe_grid(
  dim = length(hypers),
  level = c(1, 2),
  cut_off = c(0, 1.9),
  sd = sqrt(diag(sd_out$cov.fixed))
)

#' Total number of grid points is 2^1 = 2
#' Going for the minimum possible grid :)
prod(base_grid$level)

control <- aghq::default_control_tmb()

#' Start expose fit_aghq(tmb_inputs, k = k, basegrid = sparse_grid, control = control)
tmb_input <- tmb_inputs
k <- k
basegrid <- base_grid
control <- control
inner_verbose <- FALSE
progress <- NULL
map <- NULL
DLL <- "naomi_simple"

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
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

start <- Sys.time()

quad <- aghq(
  ff = ff, k = k, transformation = transformation,
  startingvalue = startingvalue, optresults = optresults,
  basegrid = basegrid, control = control
)

end <- Sys.time()

end - start
