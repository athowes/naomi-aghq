#' A local version of aghq::normalize_logpost
#' I have added an argument whereby the basegrid provided can not be adapted (if it is adapted already)
local_normalize_logpost <- function(optresults, k, whichfirst = 1, basegrid = NULL, adapt = TRUE, ndConstruction = "product", ...) {
  S <- length(optresults$mode)
  thegrid <- basegrid
  idxorder <- c(whichfirst, (1:S)[-whichfirst])

  if(adapt) {
    m <- optresults$mode[idxorder]
    H <- optresults$hessian[idxorder, idxorder]
    mvQuad::rescale(thegrid, m = m, C = Matrix::forceSymmetric(solve(H)), dec.type = 2)
  }

  nodesandweights <- cbind(mvQuad::getNodes(thegrid), mvQuad::getWeights(thegrid))
  colnames(nodesandweights) <- c(paste0("theta", idxorder), "weights")
  nodesandweights <- as.data.frame(nodesandweights)
  thetaorder <- paste0("theta", 1:S)
  if (length(idxorder) == 1) {
    nodesandweights$logpost <- sapply(nodesandweights[, thetaorder], optresults$ff$fn, ...)
  }
  else {
    nodesandweights$logpost <- apply(nodesandweights[, thetaorder], 1, optresults$ff$fn, ...)
  }
  ww <- nodesandweights$weights
  pp <- nodesandweights$logpost
  lognormconst <- aghq:::logsumexpweights(pp, ww)
  nodesandweights$logpost_normalized <- nodesandweights$logpost - lognormconst
  list(nodesandweights = nodesandweights, grid = thegrid, lognormconst = lognormconst)
}

#' A local version of aghq::aghq
#' Compatible with local_normalize_logpost -- i.e. allowing no adaption
local_aghq <- function(ff, k, startingvalue, transformation = default_transformation(), optresults = NULL, basegrid = NULL, adapt = TRUE, control = default_control(), ...) {

  validate_control(control)
  validate_transformation(transformation)
  transformation <- make_transformation(transformation)

  # If they provided a basegrid, get the k from that. If they also provided a k, compare them and issue a warning
  if (!is.null(basegrid)) {
    if (missing(k)) {
      k <- max(as.numeric(basegrid$level))
    } else {
      k2 <- max(as.numeric(basegrid$level))
      if (k != k2) {
        warning(paste0("You provided a basegrid and a specified number of quadrature points k. You do not need to specify k if you supply a basegrid. Further, they don't match: your grid has k = ",k2,", but you specified k = ",k,". Proceeding with k = ",k2,", from the supplied grid.\n"))
        k <- k2
      }
    }
  }

  # Optimization
  if (is.null(optresults)) utils::capture.output(optresults <- optimize_theta(ff, startingvalue, control, ...))

  # Normalization
  normalized_posterior <- local_normalize_logpost(optresults, k, basegrid = basegrid, adapt = adapt, ndConstruction = control$ndConstruction, ...)

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
  marginals <- vector(mode = "list", length = d)

  if (control$method_summaries[1] == "correct") {
    for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "correct")
  } else {
    for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "reuse")
  }

  out$marginals <- marginals
  out
}
