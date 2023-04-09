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

#' Debug the above
# function(ff,k,startingvalue,transformation = default_transformation(),optresults = NULL,basegrid = NULL,control = default_control(),...) {
ff <- ff
k <- k
transformation <- transformation
startingvalue <- startingvalue
optresults <- optresults
basegrid <- basegrid
control <- control

validate_control(control)
validate_transformation(transformation)
transformation <- make_transformation(transformation)
# ff <- make_moment_function(ff)
# validate_moment(ff)

print(paste0("The value of k is ", k))

# Optimization
if (is.null(optresults)) utils::capture.output(optresults <- aghq::optimize_theta(ff,startingvalue,control))

print("Optimisation complete")

# Normalization
normalized_posterior <- aghq::normalize_logpost(optresults, k, basegrid = basegrid, ndConstruction = control$ndConstruction, dec.type = dec.type)

print("Normalisation complete")

if (control$onlynormconst) return(normalized_posterior$lognormconst)

out <- list(
  normalized_posterior = normalized_posterior,
  # marginals = marginals,
  optresults = optresults,
  control = control,
  transformation = transformation
)

class(out) <- "aghq"

# Marginals

d <- length(startingvalue)
marginals <- vector(mode = "list",length = d)

if (control$method_summaries[1] == 'correct') {
  # Getting k = 2 again here
  # It's because get_... uses levels[1]
  # Change this to min(levels)
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out,j,method = 'correct')
} else {
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out,j,method = 'reuse')
}

out$marginals <- marginals
quad <- out

## Add on the info needed for marginal Laplace ##
# Add on the quad point modes and curvatures
distinctthetas <- quad$normalized_posterior$nodesandweights[ ,grep('theta',colnames(quad$normalized_posterior$nodesandweights))]
if (!is.data.frame(distinctthetas)) distinctthetas <- data.frame(theta1 = distinctthetas)

modesandhessians <- distinctthetas
if (is.null(thetanames)) {
  thetanames <- colnames(distinctthetas)
} else {
  colnames(modesandhessians) <- thetanames
  colnames(quad$normalized_posterior$nodesandweights)[grep('theta',colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
}

modesandhessians$mode <- vector(mode = 'list',length = nrow(distinctthetas))
modesandhessians$H <- vector(mode = 'list',length = nrow(distinctthetas))

for (i in 1:nrow(distinctthetas)) {
  # Get the theta
  theta <- as.numeric(modesandhessians[i,thetanames])
  # Set up the mode and hessian of the random effects. This happens when you run
  # the TMB objective with a particular theta
  ff$fn(theta)
  # Now pull the mode and hessian. Have to be careful about scoping
  mm <- ff$env$last.par
  modesandhessians[i,'mode'] <- list(list(mm[ff$env$random]))
  H <- ff$env$spHess(mm,random = TRUE)
  H <- rlang::duplicate(H) # Somehow, TMB puts all evals of spHess in the same memory location.
  modesandhessians[i,'H'] <- list(list(H))
}

quad$modesandhessians <- modesandhessians

class(quad) <- c("marginallaplace","aghq")
quad
