#' What is the delta method? What does it want? What might it know?
#' These are mysteries that have puzzled mankind (me) for too long. No more

TMB::sdreport <- function(
    obj,
    par.fixed = NULL,
    hessian.fixed = NULL,
    getJointPrecision = FALSE,
    bias.correct = FALSE,
    bias.correct.control = list(sd = FALSE, split = NULL, nsplit = NULL),
    ignore.parm.uncertainty = FALSE,
    getReportCovariance = TRUE,
    skip.delta.method = FALSE
) {

if (is.null(obj$env$ADGrad) & (!is.null(obj$env$random)))
  stop("Cannot calculate sd's without type ADGrad available in object for random effect models.")
obj2 <- MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun", ADreport = TRUE, DLL = obj$env$DLL, silent = obj$env$silent)
r <- obj$env$random
if (is.null(par.fixed)) {
  par <- obj$env$last.par.best
  if (!is.null(r))
    par.fixed <- par[-r]
  else par.fixed <- par
  gradient.fixed <- obj$gr(par.fixed)
}
else {
  gradient.fixed <- obj$gr(par.fixed)
  par <- obj$env$last.par
}
if (length(par.fixed) == 0)
  ignore.parm.uncertainty <- TRUE
if (ignore.parm.uncertainty) {
  hessian.fixed <- NULL
  pdHess <- TRUE
  Vtheta <- matrix(0, length(par.fixed), length(par.fixed))
}
else {
  if (is.null(hessian.fixed)) {
    hessian.fixed <- optimHess(par.fixed, obj$fn, obj$gr)
  }
  pdHess <- !is.character(try(chol(hessian.fixed), silent = TRUE))
  Vtheta <- try(solve(hessian.fixed), silent = TRUE)
  if (is(Vtheta, "try-error"))
    Vtheta <- hessian.fixed * NaN
}
if (!is.null(r)) {
  hessian.random <- obj$env$spHess(par, random = TRUE)
  L <- obj$env$L.created.by.newton
  if (!is.null(L)) {
    updateCholesky(L, hessian.random)
    hessian.random@factors <- list(SPdCholesky = L)
  }
}
phi <- try(obj2$fn(par), silent = TRUE)
if (is.character(phi) | length(phi) == 0) {
  phi <- numeric(0)
}
ADGradForward0Initialized <- FALSE
ADGradForward0Initialize <- function() {
  obj$env$f(par, order = 0, type = "ADGrad")
  ADGradForward0Initialized <<- TRUE
}
doDeltaMethod <- function(chunk = NULL) {
  simpleCase <- is.null(r)
  if (length(phi) == 0) {
    simpleCase <- TRUE
  }
  else {
    if (is.null(chunk)) {
      Dphi <- obj2$gr(par)
    }
    else {
      w <- rep(0, length(phi))
      phiDeriv <- function(i) {
        w[i] <- 1
        obj2$env$f(par, order = 1, rangeweight = w,
                   doforward = 0)
      }
      Dphi <- t(sapply(chunk, phiDeriv))
      phi <- phi[chunk]
    }
    if (!is.null(r)) {
      Dphi.random <- Dphi[, r, drop = FALSE]
      Dphi.fixed <- Dphi[, -r, drop = FALSE]
      if (all(Dphi.random == 0)) {
        simpleCase <- TRUE
        Dphi <- Dphi.fixed
      }
    }
  }
  if (simpleCase) {
    if (length(phi) > 0) {
      cov <- Dphi %*% Vtheta %*% t(Dphi)
    }
    else cov <- matrix(, 0, 0)
  }
  else {
    tmp <- solve(hessian.random, t(Dphi.random))
    tmp <- as.matrix(tmp)
    term1 <- Dphi.random %*% tmp
    if (ignore.parm.uncertainty) {
      term2 <- 0
    }
    else {
      f <- obj$env$f
      w <- rep(0, length(par))
      if (!ADGradForward0Initialized)
        ADGradForward0Initialize()
      reverse.sweep <- function(i) {
        w[r] <- tmp[, i]
        -f(par, order = 1, type = "ADGrad", rangeweight = w,
           doforward = 0)[-r]
      }
      A <- t(do.call("cbind", lapply(seq_along(phi),
                                     reverse.sweep))) + Dphi.fixed
      term2 <- A %*% (Vtheta %*% t(A))
    }
    cov <- term1 + term2
  }
  cov
}
if (!skip.delta.method) {
  if (getReportCovariance) {
    cov <- doDeltaMethod()
    sd <- sqrt(diag(cov))
  }
  else {
    tmp <- lapply(seq_along(phi), doDeltaMethod)
    sd <- sqrt(as.numeric(unlist(tmp)))
    cov <- NA
  }
}
else {
  sd <- rep(NA, length(phi))
  cov <- NA
}
ans <- list(value = phi, sd = sd, cov = cov, par.fixed = par.fixed,
            cov.fixed = Vtheta, pdHess = pdHess, gradient.fixed = gradient.fixed)
if (bias.correct) {
  epsilon <- rep(0, length(phi))
  names(epsilon) <- names(phi)
  parameters <- obj$env$parameters
  parameters$TMB_epsilon_ <- epsilon
  doEpsilonMethod <- function(chunk = NULL) {
    if (!is.null(chunk)) {
      mapfac <- rep(NA, length(phi))
      mapfac[chunk] <- chunk
      parameters$TMB_epsilon_ <- updateMap(parameters$TMB_epsilon_,
                                           factor(mapfac))
    }
    obj3 <- MakeADFun(obj$env$data, parameters, random = obj$env$random,
                      checkParameterOrder = FALSE, DLL = obj$env$DLL,
                      silent = obj$env$silent)
    obj3$env$start <- c(par, epsilon)
    obj3$env$random.start <- expression(start[random])
    h <- obj$env$spHess(random = TRUE)
    h3 <- obj3$env$spHess(random = TRUE)
    pattern.unchanged <- identical(h@i, h3@i) & identical(h@p,
                                                          h3@p)
    if (pattern.unchanged) {
      if (!obj$env$silent)
        cat("Re-using symbolic Cholesky\n")
      obj3$env$L.created.by.newton <- L
    }
    else {
      if (.Call("have_tmb_symbolic", PACKAGE = "TMB"))
        runSymbolicAnalysis(obj3)
    }
    if (!is.null(chunk))
      epsilon <- epsilon[chunk]
    par.full <- c(par.fixed, epsilon)
    i <- (1:length(par.full)) > length(par.fixed)
    grad <- obj3$gr(par.full)
    Vestimate <- if (bias.correct.control$sd) {
      hess <- numDeriv::jacobian(obj3$gr, par.full)
      -hess[i, i] + hess[i, !i] %*% Vtheta %*% hess[!i,
                                                    i]
    }
    else matrix(NA)
    estimate <- grad[i]
    names(estimate) <- names(epsilon)
    list(value = estimate, sd = sqrt(diag(Vestimate)),
         cov = Vestimate)
  }
  nsplit <- bias.correct.control$nsplit
  if (is.null(nsplit)) {
    split <- bias.correct.control$split
  }
  else {
    split <- split(seq_along(phi), cut(seq_along(phi),
                                       nsplit))
  }
  if (is.null(split)) {
    ans$unbiased <- doEpsilonMethod()
  }
  else {
    tmp <- lapply(split, doEpsilonMethod)
    m <- if (bias.correct.control$sd)
      length(phi)
    else 1
    ans$unbiased <- list(value = rep(NA, length(phi)),
                         sd = rep(NA, m), cov = matrix(NA, m, m))
    for (i in seq_along(split)) {
      ans$unbiased$value[split[[i]]] <- tmp[[i]]$value
      if (bias.correct.control$sd) {
        ans$unbiased$sd[split[[i]]] <- tmp[[i]]$sd
        ans$unbiased$cov[split[[i]], split[[i]]] <- tmp[[i]]$cov
      }
    }
  }
}
if (!is.null(r)) {
  if (is(L, "dCHMsuper")) {
    diag.term1 <- solveSubset(L = L, diag = TRUE)
    if (ignore.parm.uncertainty) {
      diag.term2 <- 0
    }
    else {
      f <- obj$env$f
      w <- rep(0, length(par))
      if (!ADGradForward0Initialized)
        ADGradForward0Initialize()
      reverse.sweep <- function(i) {
        w[i] <- 1
        f(par, order = 1, type = "ADGrad", rangeweight = w,
          doforward = 0)[r]
      }
      nonr <- setdiff(seq_along(par), r)
      framework <- .Call("getFramework", PACKAGE = obj$env$DLL)
      if (framework != "TMBad")
        tmp <- sapply(nonr, reverse.sweep)
      else tmp <- f(par, order = 1, type = "ADGrad",
                    keepx = nonr, keepy = r)
      if (!is.matrix(tmp))
        tmp <- matrix(tmp, ncol = length(nonr))
      A <- solve(hessian.random, tmp)
      diag.term2 <- rowSums((A %*% Vtheta) * A)
    }
    ans$par.random <- par[r]
    ans$diag.cov.random <- diag.term1 + diag.term2
    if (getJointPrecision) {
      if (length(par.fixed) == 0) {
        ans$jointPrecision <- hessian.random
      }
      else if (!ignore.parm.uncertainty) {
        G <- hessian.random %*% A
        G <- as.matrix(G)
        M1 <- cbind2(hessian.random, G)
        M2 <- cbind2(t(G), as.matrix(t(A) %*% G) +
                       hessian.fixed)
        M <- rbind2(M1, M2)
        M <- forceSymmetric(M, uplo = "L")
        dn <- c(names(par)[r], names(par[-r]))
        dimnames(M) <- list(dn, dn)
        p <- invPerm(c(r, (1:length(par))[-r]))
        ans$jointPrecision <- M[p, p]
      }
      else {
        warning("ignore.parm.uncertainty ==> No joint precision available")
      }
    }
  }
  else {
    warning("Could not report sd's of full randomeffect vector.")
  }
}
ans$env <- new.env(parent = emptyenv())
ans$env$parameters <- obj$env$parameters
ans$env$random <- obj$env$random
ans$env$ADreportDims <- obj2$env$ADreportDims
class(ans) <- "sdreport"
ans
}
