# sparse_quad <- fit_aghq(tmb_inputs, k = 2, basegrid = sparse_grid, control = control)
# fit_aghq <- function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, ...)

tmb_input <- tmb_inputs
inner_verbose <- FALSE
progress <- NULL
map <- NULL

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map)

# quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = 2, basegrid = sparse_grid, control = control)

ff <- obj
k <- 2
startingvalue <- obj$par
transformation <- default_transformation()
optresults <- NULL
basegrid <- sparse_grid
control <- control

validate_control(control, type = "tmb")
validate_transformation(transformation)
transformation <- make_transformation(transformation)
thetanames <- NULL
if (exists("par", ff))
  thetanames <- make.unique(names(ff$par), sep = "")
if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr, theta, method = "Richardson")
}

# quad <- aghq(ff = ff, k = k, transformation = transformation, startingvalue = startingvalue, optresults = optresults, basegrid = basegrid, control = control)

validate_control(control)
validate_transformation(transformation)
transformation <- make_transformation(transformation)

if (!is.null(basegrid)) {
  if (missing(k)) {
    k <- max(as.numeric(basegrid$level))
  }
  else {
    k2 <- max(as.numeric(basegrid$level))
    if (k != k2) {
      warning(paste0("You provided a basegrid and a specified number of quadrature points k. You do not need to specify k if you supply a basegrid. Further, they don't match: your grid has k = ",
                     k2, ", but you specified k = ", k, ". Proceeding with k = ",
                     k2, ", from the supplied grid.\n"))
      k <- k2
    }
  }
}

if (is.null(optresults)) utils::capture.output(optresults <- optimize_theta(ff, startingvalue, control))

normalized_posterior <- normalize_logpost(optresults, k, basegrid = basegrid, ndConstruction = control$ndConstruction)

if (control$onlynormconst) return(normalized_posterior$lognormconst)

out <- list(normalized_posterior = normalized_posterior, optresults = optresults, control = control, transformation = transformation)

class(out) <- "aghq"
d <- length(startingvalue)
marginals <- vector(mode = "list", length = d)

if (control$method_summaries[1] == "correct") {
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "correct")
} else {
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "reuse")
}

out$marginals <- marginals
out

method <- method[1]
if (method == "auto") method <- "reuse"
if (method == "reuse") return(marginal_posterior.list(quad$optresults, get_numquadpoints(quad), j, ...))

S <- get_param_dim(quad)
idxorder <- c(j, (1:S)[-j])
thetaminusj <- (1:S)[-j]
if (is.null(qq)) {
  mm <- get_mode(quad)[idxorder]
  HH <- get_hessian(quad)[idxorder, idxorder]
  gg <- mvQuad::createNIGrid(1, "GHe", get_numquadpoints(quad))
  mvQuad::rescale(gg, m = mm[1], C = solve(HH)[1, 1], dec.type = 2)
  qqq <- as.numeric(mvQuad::getNodes(gg))
  out <- vector(mode = "list", length = length(qqq))
  for (i in 1:length(qqq)) {
    out[[i]] <- marginal_posterior.aghq(quad = quad, j = j, qq = qqq[i], method = "correct", ...)
  }
  out <- Reduce(rbind, out)
} else {
  cname <- colnames(get_nodesandweights(quad))[j]
  if (S == 1) {
    out <- data.frame(qq, quad$optresults$ff$fn(qq) -
                        get_log_normconst(quad))
    colnames(out) <- c(cname, "logmargpost")
    return(out)
  }
  fn <- function(theta) quad$optresults$ff$fn(splice(theta,
                                                     qq, j))
  gr <- function(theta) quad$optresults$ff$gr(splice(theta,
                                                     qq, j))[-j]
  he <- function(theta) quad$optresults$ff$he(splice(theta,
                                                     qq, j))[-j, -j]
  ffm <- list(fn = fn, gr = gr, he = he)
  newcontrol <- quad$control
  newcontrol$onlynormconst <- TRUE
  newcontrol$negate <- FALSE
  lognumerator <- get_log_normconst(aghq(ffm, get_numquadpoints(quad),
                                         get_mode(quad)[-j], control = newcontrol))
  out <- data.frame(qq, lognumerator - get_log_normconst(quad))
  colnames(out) <- c(cname, "logmargpost")
}
out

if (control$onlynormconst) return(quad)

distinctthetas <- quad$normalized_posterior$nodesandweights[, grep("theta", colnames(quad$normalized_posterior$nodesandweights))]

if (!is.data.frame(distinctthetas))
  distinctthetas <- data.frame(theta1 = distinctthetas)
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

quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")

quad$obj <- obj
quad
