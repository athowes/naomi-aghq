tmb_input <- tmb_inputs
k <- 2
basegrid <- sparse_grid

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress)

ff <- obj
optresults <- NULL
startingvalue <- obj$par

transformation = default_transformation()
control = default_control_tmb()

validate_control(control, type = "tmb")
validate_transformation(transformation)
transformation <- make_transformation(transformation)
thetanames <- NULL
if (exists("par", ff))
  thetanames <- make.unique(names(ff$par), sep = "")
if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr, theta,
                                              method = "Richardson")
}

quad <- aghq(ff = ff, k = k, transformation = transformation,
             startingvalue = startingvalue, optresults = optresults,
             basegrid = basegrid, control = control)

#' It's this quad which crashes! Expose this

if (control$onlynormconst)
  return(quad)
distinctthetas <- quad$normalized_posterior$nodesandweights[,
                                                            grep("theta", colnames(quad$normalized_posterior$nodesandweights))]
if (!is.data.frame(distinctthetas))
  distinctthetas <- data.frame(theta1 = distinctthetas)
modesandhessians <- distinctthetas
if (is.null(thetanames)) {
  thetanames <- colnames(distinctthetas)
}
else {
  colnames(modesandhessians) <- thetanames
  colnames(quad$normalized_posterior$nodesandweights)[grep("theta",
                                                           colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
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
quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")
quad
