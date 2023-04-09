k <- 2
(levels <- c(rep(k, s), rep(1, n_hyper - s)))
(pca_base_grid <- mvQuad::createNIGrid(dim = n_hyper, type = "GHe", level = levels))

#' Test fitting without writing basegrid: runs fine
quad <- fit_aghq(tmb_inputs_simple, k = 1, dec.type = 1)
stopifnot(nrow(quad$modesandhessians) == 1)

#' Test fitting with writing basegrid
quad <- fit_aghq(tmb_inputs_simple, k = 1, basegrid = pca_base_grid, dec.type = 1)

#' Debug the above
# function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, DLL = "naomi_simple", ...) {
tmb_input <- tmb_inputs_simple
inner_verbose <- FALSE
progress <- NULL
map <- NULL
DLL <- "naomi_simple"
k <- 1
basegrid <- pca_base_grid
dec.type <- 1

if (DLL == "naomi_simple") {
  stopifnot(inherits(tmb_input, "naomi_simple_tmb_input"))
} else {
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
}

obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL = DLL)

quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = k, basegrid = basegrid, dec.type = dec.type)

#' Debug the above
# function(ff, k, startingvalue, transformation = default_transformation(), optresults = NULL, basegrid = NULL, control = default_control_tmb(), ...) {
ff <- obj
k <- k
startingvalue <- obj$par
basegrid <- basegrid
dec.type <- dec.type
transformation <- default_transformation()
control <- default_control_tmb()
optresults <- NULL

validate_control(control,type='tmb')
validate_transformation(transformation)
transformation <- make_transformation(transformation)

thetanames <- NULL
if (exists('par',ff)) thetanames <- make.unique(names(ff$par),sep='')
if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr,theta,method = 'Richardson')
}

#' This runs
# quad <- aghq(ff = ff, k = k, transformation = transformation, startingvalue = startingvalue, optresults = optresults, basegrid = NULL, control = control)

#' This doesn't
quad <- aghq(ff = ff, k = k, transformation = transformation, startingvalue = startingvalue, optresults = optresults, basegrid = basegrid, control = control)
