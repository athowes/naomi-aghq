#' Expose internals of `MakeADFun` to look to change the matrix algebra
#' This would be a next step to implementing Wood SINLA Laplace marginals

local_MakeADFun <- function(
  data,
  parameters,
  map = list(),
  type = c("ADFun", "Fun", "ADGrad"[!intern && (!is.null(random) || !is.null(profile))]),
  random = NULL,
  profile = NULL,
  random.start = expression(last.par.best[random]),
  hessian = FALSE,
  method = "BFGS",
  inner.method = "newton",
  inner.control = list(maxit = 1000),
  MCcontrol = list(doMC = FALSE, seed = 123, n = 100),
  ADreport = FALSE,
  atomic = TRUE,
  LaplaceNonZeroGradient = FALSE,
  DLL = getUserDLL(),
  checkParameterOrder = TRUE,
  regexp = FALSE,
  silent = FALSE,
  intern = FALSE,
  integrate = NULL,
  ...) {

  if (!DLL %in% names(getLoadedDLLs())) {
    stop(sprintf("'%s' was not found in the list of loaded DLLs. Forgot to dyn.load(dynlib('%s')) ?", DLL, DLL))
  }

  env <- environment()

  if (!is.list(data))
    stop("'data' must be a list")

  ok <- function(x) (is.matrix(x) | is.vector(x) | is.array(x)) & (is.numeric(x) | is.logical(x))

  ok.data <- function(x) ok(x) | is.factor(x) | is(x, "sparseMatrix") | is.list(x) | (is.character(x) & length(x) == 1)

  check.passed <- function(x) {
    y <- attr(x, "check.passed")
    if (is.null(y)) FALSE
    else y
  }

  if (!check.passed(data)) {
    if (!all(sapply(data, ok.data))) {
      cat("Problem with these data entries:\n")
      print(which(!sapply(data, ok.data)))
      stop("Only numeric matrices, vectors, arrays, ",
           "factors, lists or length-1-characters ",
           "can be interfaced")
    }
  }

  if (!check.passed(parameters)) {
    if (!all(sapply(parameters, ok))) {
      cat("Problem with these parameter entries:\n")
      print(which(!sapply(parameters, ok)))
      stop("Only numeric matrices, vectors and arrays ",
           "can be interfaced")
    }
  }

  if (length(data)) {
    dataSanitize <- function(x) {
      if (is.list(x))
        return(lapply(x, dataSanitize))
      if (is(x, "sparseMatrix")) {
        x <- as(x, "TsparseMatrix")
      }
      else if (is.character(x)) {
      }
      else {
        if (is.factor(x))
          x <- unclass(x) - 1L
        storage.mode(x) <- "double"
      }
      x
    }
    if (!check.passed(data)) {
      data <- lapply(data, dataSanitize)
    }
    attr(data, "check.passed") <- TRUE
  }

  if (length(parameters)) {
    parameterSanitize <- function(x) {
      storage.mode(x) <- "double"
      x
    }
    if (!check.passed(parameters)) {
      parameters <- lapply(parameters, parameterSanitize)
    }
    attr(parameters, "check.passed") <- TRUE
  }

  if (checkParameterOrder) {
    parNameOrder <- getParameterOrder(data, parameters, new.env(),
                                      DLL = DLL)
    if (!identical(names(parameters), parNameOrder)) {
      if (!silent)
        cat("Order of parameters:\n")
      if (!silent)
        print(names(parameters))
      if (!silent)
        cat("Not matching template order:\n")
      if (!silent)
        print(parNameOrder)
      keepAttrib(parameters) <- parameters[parNameOrder]
      if (!silent)
        cat("Your parameter list has been re-ordered.\n(Disable this warning with checkParameterOrder=FALSE)\n")
    }
  }

  if (length(map) > 0) {
    ok <- all(names(map) %in% names(parameters))
    if (!ok)
      stop("Names in map must correspond to parameter names")
    ok <- all(sapply(map, is.factor))
    if (!ok)
      stop("map must contain factors")
    ok <- sapply(parameters[names(map)], length) == sapply(map,
                                                           length)
    if (!all(ok))
      stop("A map factor length must equal parameter length")
    param.map <- lapply(names(map), function(nam) {
      updateMap(parameters[[nam]], map[[nam]])
    })
    keepAttrib(parameters[names(map)]) <- param.map
  }
  lrandom <- function() {
    ans <- logical(length(par))
    ans[random] <- TRUE
    ans
  }
  lfixed <- function() {
    !lrandom()
  }
  parList <- function(x = par[lfixed()], par = last.par) {
    ans <- parameters
    nonemp <- sapply(ans, function(x) length(x) > 0)
    nonempindex <- which(nonemp)
    skeleton <- as.relistable(ans[nonemp])
    par[lfixed()] <- x
    li <- relist(par, skeleton)
    reshape <- function(x) {
      if (is.null(attr(x, "map")))
        return(x)
      y <- attr(x, "shape")
      f <- attr(x, "map")
      i <- which(f >= 0)
      y[i] <- x[f[i] + 1]
      y
    }
    for (i in seq(skeleton)) {
      ans[[nonempindex[i]]][] <- as.vector(li[[i]])
    }
    for (i in seq(ans)) {
      ans[[i]] <- reshape(ans[[i]])
    }
    ans
  }
  type <- match.arg(type, eval(type), several.ok = TRUE)
  reportenv <- new.env()
  par <- NULL
  last.par.ok <- last.par <- last.par1 <- last.par2 <- last.par.best <- NULL
  value.best <- Inf
  ADFun <- NULL
  Fun <- NULL
  ADGrad <- NULL
  tracepar <- FALSE
  validpar <- function(x) TRUE
  tracemgc <- TRUE
  L.created.by.newton <- skipFixedEffects <- spHess <- NULL
  beSilent <- function() {
    tracemgc <<- FALSE
    inner.control$trace <<- FALSE
    silent <<- TRUE
    cf <- config(DLL = DLL)
    i <- grep("^trace.", names(cf))
    cf[i] <- 0
    cf$DLL <- DLL
    do.call(config, cf)
    NULL
  }
  if (silent)
    beSilent()
  ADreportDims <- NULL
  ADreportIndex <- function() {
    lngt <- sapply(ADreportDims, prod)
    offset <- head(cumsum(c(1, lngt)), -1)
    ans <- lapply(seq_along(lngt), function(i) array(seq(from = offset[i],
                                                         length.out = lngt[i]), ADreportDims[[i]]))
    names(ans) <- names(ADreportDims)
    ans
  }
  .random <- random
  retape <- function(set.defaults = TRUE) {
    omp <- config(DLL = DLL)
    random <<- .random
    if (atomic) {
      Fun <<- MakeDoubleFunObject(data, parameters, reportenv,
                                  DLL = DLL)
      out <- EvalDoubleFunObject(Fun, unlist(parameters),
                                 get_reportdims = TRUE)
      ADreportDims <<- attr(out, "reportdims")
    }
    if (is.character(profile)) {
      random <<- c(random, profile)
    }
    if (is.character(random)) {
      if (!regexp) {
        if (!all(random %in% names(parameters))) {
          cat("Some 'random' effect names does not match 'parameter' list:\n")
          print(setdiff(random, names(parameters)))
          cat("(Note that regular expression match is disabled by default)\n")
          stop()
        }
        if (any(duplicated(random))) {
          cat("Duplicates in 'random' - will be removed\n")
          random <<- unique(random)
        }
        tmp <- lapply(parameters, function(x) x * 0)
        tmp[random] <- lapply(tmp[random], function(x) x *
                                0 + 1)
        random <<- which(as.logical(unlist(tmp)))
        if (length(random) == 0)
          random <<- NULL
      }
      if (regexp) {
        random <<- grepRandomParameters(parameters, random)
        if (length(random) == 0) {
          cat("Selected random effects did not match any model parameters.\n")
          random <<- NULL
        }
      }
      if (is.character(profile)) {
        tmp <- lapply(parameters, function(x) x * 0)
        tmp[profile] <- lapply(tmp[profile], function(x) x *
                                 0 + 1)
        profile <<- match(which(as.logical(unlist(tmp))),
                          random)
        if (length(profile) == 0)
          random <<- NULL
        if (any(duplicated(profile)))
          stop("Profile parameter vector not unique.")
        tmp <- rep(0L, length(random))
        tmp[profile] <- 1L
        profile <<- tmp
      }
      if (set.defaults) {
        par <<- unlist(parameters)
      }
    }
    if ("ADFun" %in% type) {
      if (omp$autopar)
        openmp(1, DLL = DLL)
      ADFun <<- MakeADFunObject(data, parameters, reportenv,
                                ADreport = ADreport, DLL = DLL)
      if (omp$autopar)
        openmp(omp$nthreads, DLL = DLL)
      if (!is.null(integrate)) {
        nm <- sapply(parameters, length)
        nmpar <- rep(names(nm), nm)
        for (i in seq_along(integrate)) {
          I <- integrate[i]
          if (is.null(names(I)) || names(I) == "") {
            I <- I[[1]]
          }
          ok <- all(names(I) %in% nmpar[random])
          if (!ok)
            stop("Names to be 'integrate'd must be among the random parameters")
          w <- which(nmpar[random] %in% names(I))
          arg_which <- I[[1]]$which
          if (!is.null(arg_which))
            w <- w[arg_which]
          method <- sapply(I, function(x) x$method)
          ok <- all(duplicated(method)[-1])
          if (!ok)
            stop("Grouping only allowed for identical methods")
          method <- method[1]
          cfg <- NULL
          if (method == "marginal_sr") {
            fac <- factor(nmpar[random[w]], levels = names(I))
            cfg <- list(grid = I, random2grid = fac)
          }
          else {
            cfg <- I[[1]]
          }
          stopifnot(is.list(cfg))
          TransformADFunObject(ADFun, method = method,
                               random_order = random[w], config = cfg, mustWork = 1L)
          activeDomain <- as.logical(info(ADFun)$activeDomain)
          random_remove <- random[w][!activeDomain[random[w]]]
          TransformADFunObject(ADFun, method = "remove_random_parameters",
                               random_order = random_remove, mustWork = 1L)
          attr(ADFun$ptr, "par") <- attr(ADFun$ptr, "par")[-random_remove]
          par_mask <- rep(FALSE, length(attr(ADFun$ptr,
                                             "par")))
          par_mask[random] <- TRUE
          par <<- par[-random_remove]
          nmpar <- nmpar[-random_remove]
          par_mask <- par_mask[-random_remove]
          random <<- which(par_mask)
          if (length(random) == 0) {
            random <<- NULL
            type <<- setdiff(type, "ADGrad")
          }
          if (config(DLL = DLL)$optimize.instantly) {
            TransformADFunObject(ADFun, method = "optimize",
                                 mustWork = 1L)
          }
        }
      }
      if (intern) {
        cfg <- inner.control
        if (is.null(cfg$sparse))
          cfg$sparse <- TRUE
        cfg <- lapply(cfg, as.double)
        TransformADFunObject(ADFun, method = "laplace",
                             config = cfg, random_order = random, mustWork = 1L)
        TransformADFunObject(ADFun, method = "remove_random_parameters",
                             random_order = random, mustWork = 1L)
        attr(ADFun$ptr, "par") <- attr(ADFun$ptr, "par")[-random]
        par <<- par[-random]
        random <<- NULL
        if (config(DLL = DLL)$optimize.instantly) {
          TransformADFunObject(ADFun, method = "optimize",
                               mustWork = 1L)
        }
      }
      if (set.defaults) {
        par <<- attr(ADFun$ptr, "par")
        last.par <<- par
        last.par1 <<- par
        last.par2 <<- par
        last.par.best <<- par
        value.best <<- Inf
      }
    }
    if (omp$autopar && !ADreport) {
      TransformADFunObject(ADFun, method = "parallel_accumulate",
                           num_threads = as.integer(openmp(DLL = DLL)),
                           mustWork = 0L)
    }
    if (length(random) > 0) {
      TransformADFunObject(ADFun, method = "reorder_random",
                           random_order = random, mustWork = 0L)
    }
    if ("Fun" %in% type) {
      Fun <<- MakeDoubleFunObject(data, parameters, reportenv,
                                  DLL = DLL)
    }
    if ("ADGrad" %in% type) {
      retape_adgrad(lazy = TRUE)
    }
    env$skipFixedEffects <- !is.null(ADGrad)
    delayedAssign("spHess", sparseHessianFun(env, skipFixedEffects = skipFixedEffects),
                  assign.env = env)
  }
  retape_adgrad <- function(lazy = TRUE) {
    if (!lazy)
      random <- NULL
    ADGrad <<- MakeADGradObject(data, parameters, reportenv,
                                random = random, f = ADFun$ptr, DLL = DLL)
  }
  retape(set.defaults = TRUE)
  usingAtomics <- function() .Call("usingAtomics", PACKAGE = DLL)
  .data <- NULL
  f <- function(theta = par, order = 0, type = "ADdouble",
                cols = NULL, rows = NULL, sparsitypattern = 0, rangecomponent = 1,
                rangeweight = NULL, dumpstack = 0, doforward = 1, do_simulate = 0,
                set_tail = 0, keepx = NULL, keepy = NULL) {
    if (isNullPointer(ADFun$ptr)) {
      if (silent)
        beSilent()
      retape(set.defaults = FALSE)
    }
    data_changed <- !identical(.data, data)
    if (data_changed) {
      .data <<- data
    }
    switch(type, ADdouble = {
      res <- EvalADFunObject(ADFun, theta, order = order,
                             hessiancols = cols, hessianrows = rows, sparsitypattern = sparsitypattern,
                             rangecomponent = rangecomponent, rangeweight = rangeweight,
                             dumpstack = dumpstack, doforward = doforward,
                             set_tail = set_tail, data_changed = data_changed)
      last.par <<- theta
      if (order == 1) last.par1 <<- theta
      if (order == 2) last.par2 <<- theta
    }, double = {
      res <- EvalDoubleFunObject(Fun, theta, do_simulate = do_simulate)
    }, ADGrad = {
      res <- EvalADFunObject(ADGrad, theta, order = order,
                             hessiancols = cols, hessianrows = rows, sparsitypattern = sparsitypattern,
                             rangecomponent = rangecomponent, rangeweight = rangeweight,
                             dumpstack = dumpstack, doforward = doforward,
                             set_tail = set_tail, keepx = keepx, keepy = keepy,
                             data_changed = data_changed)
    }, stop("invalid 'type'"))
    res
  }
  h <- function(theta = par, order = 0, hessian, L, ...) {
    if (order == 0) {
      logdetH <- 2 * determinant(L)$mod
      ans <- f(theta, order = 0) + 0.5 * logdetH - length(random)/2 *
        log(2 * pi)
      if (LaplaceNonZeroGradient) {
        grad <- f(theta, order = 1)[random]
        ans - 0.5 * sum(grad * as.numeric(solveCholesky(L,
                                                        grad)))
      }
      else ans
    }
    else if (order == 1) {
      if (LaplaceNonZeroGradient)
        stop("Not correct for LaplaceNonZeroGradient=TRUE")
      e <- environment(spHess)
      solveSubset <- function(L) .Call("tmb_invQ", L, PACKAGE = "TMB")
      solveSubset2 <- function(L) .Call("tmb_invQ_tril_halfdiag",
                                        L, PACKAGE = "TMB")
      ihessian <- solveSubset2(L)
      if (!is.null(profile)) {
        perm <- L@perm + 1L
        ihessian <- .Call("tmb_sparse_izamd", ihessian,
                          profile[perm], 0, PACKAGE = "TMB")
      }
      lookup <- function(A, B, r = NULL) {
        A <- tril(A)
        B <- tril(B)
        B@x[] <- seq.int(length.out = length(B@x))
        if (!is.null(r)) {
          B <- .Call("tmb_half_diag", B, PACKAGE = "TMB")
          B <- tril(B[r, r, drop = FALSE]) + tril(t(B)[r,
                                                       r, drop = FALSE])
        }
        m <- .Call("match_pattern", A, B, PACKAGE = "TMB")
        B@x[m]
      }
      if (is.null(e$ind1)) {
        if (!silent)
          cat("Matching hessian patterns... ")
        iperm <- invPerm(L@perm + 1L)
        e$ind1 <- lookup(hessian, ihessian, iperm)
        e$ind2 <- lookup(hessian, e$Hfull, random)
        if (!silent)
          cat("Done\n")
      }
      w <- rep(0, length = length(e$Hfull@x))
      w[e$ind2] <- ihessian@x[e$ind1]
      as.vector(f(theta, order = 1)) + EvalADFunObject(e$ADHess,
                                                       theta, order = 1, rangeweight = w)
    }
    else stop(sprintf("'order'=%d not yet implemented", order))
  }
  ff <- function(par.fixed = par[-random], order = 0, ...) {
    names(par.fixed) <- names(par[-random])
    f0 <- function(par.random, order = 0, ...) {
      par[random] <- par.random
      par[-random] <- par.fixed
      res <- f(par, order = order, set_tail = random[1],
               ...)
      switch(order + 1, res, res[random], res[random, random])
    }
    H0 <- function(par.random) {
      par[random] <- par.random
      par[-random] <- par.fixed
      spHess(par, random = TRUE, set_tail = random[1])
    }
    if (inner.method == "newton") {
      opt <- try(do.call("newton", c(list(par = eval(random.start),
                                          fn = f0, gr = function(x) f0(x, order = 1), he = H0,
                                          env = env), inner.control)), silent = silent)
      if (inherits(opt, "try-error") || !is.finite(opt$value)) {
        if (order == 0)
          return(NaN)
        if (order == 1)
          stop("inner newton optimization failed during gradient calculation")
        stop("invalid 'order'")
      }
    }
    else {
      opt <- optim(eval(random.start), fn = f0, gr = function(x) f0(x,
                                                                    order = 1), method = inner.method, control = inner.control)
    }
    par[random] <- opt$par
    par[-random] <- par.fixed
    if (!skipFixedEffects) {
      hess <- spHess(par)
      hessian <- hess[random, random]
    }
    else {
      hessian <- spHess(par, random = TRUE)
    }
    if (!is.null(profile)) {
      hessian <- .Call("tmb_sparse_izamd", hessian, profile,
                       1, PACKAGE = "TMB")
    }
    if (inherits(env$L.created.by.newton, "dCHMsuper")) {
      L <- env$L.created.by.newton
      updateCholesky(L, hessian)
    }
    else L <- Cholesky(hessian, perm = TRUE, LDL = FALSE,
                       super = TRUE)
    if (order == 0) {
      res <- h(par, order = 0, hessian = hessian, L = L)
      if (!is.null(profile)) {
        res <- res + sum(profile)/2 * log(2 * pi)
      }
      if (is.finite(res)) {
        if (res < value.best) {
          last.par.best <<- par
          value.best <<- res
        }
      }
    }
    if (order == 1) {
      grad <- h(par, order = 1, hessian = hessian, L = L)
      if (!is.null(profile)) {
        if (!skipFixedEffects) {
          hess <- spHess(par)
          hessian <- hess[random, random]
        }
        else {
          hessian <- spHess(par, random = TRUE)
        }
        updateCholesky(L, hessian)
      }
      if (!skipFixedEffects) {
        res <- grad[-random] - hess[-random, random] %*%
          as.vector(solveCholesky(L, grad[random]))
      }
      else {
        w <- rep(0, length(par))
        w[random] <- as.vector(solveCholesky(L, grad[random]))
        res <- grad[-random] - f(par, order = 1, type = "ADGrad",
                                 rangeweight = w)[-random]
      }
      res <- drop(res)
    }
    if (order == 2) {
      n <- length(par)
      nr <- length(random)
      nf <- n - nr
      fixed <- setdiff(1:n, random)
      D1h <- h(par, order = 1)
      D2h <- h(par, order = 2)
      D2f <- f(par, order = 2, cols = random)
      D3f <- sapply(random, function(i) f(par, type = "ADGrad",
                                          order = 2, rangecomponent = i))
      I.D2f <- solve(D2f[random, ])
      D1eta <- -t(D2f[-random, ] %*% I.D2f)
      D3f.D1eta <- D3f %*% D1eta
      dim(D3f.D1eta) <- c(n, n, nf)
      dim(D3f) <- c(n, n, nr)
      D3f.fixed <- D3f[fixed, , ]
      D2eta <- sapply(1:nf, function(i) {
        -I.D2f %*% (t(D3f.fixed[i, fixed, ]) + D3f.D1eta[random,
                                                         fixed, i] + (D3f.fixed[i, random, ] + D3f.D1eta[random,
                                                                                                         random, i]) %*% D1eta)
      })
      dim(D2eta) <- c(nr, nf, nf)
      D2h.fixed <- D2h[fixed, ]
      res <- sapply(1:nf, function(i) {
        D2h.fixed[i, fixed] + t(D2h.fixed[, random] %*%
                                  D1eta[, i]) + (t(D2h.fixed[i, random]) + t(D2h[random,
                                                                                 random] %*% D1eta[, i])) %*% D1eta + D1h[,
                                                                                                                          random] %*% D2eta[, , i]
      })
      attr(res, "D2eta") <- D2eta
      attr(res, "D1eta") <- D1eta
    }
    if (all(is.finite(res)))
      last.par.ok <<- par
    res
  }
  MC <- function(par = last.par, par0 = last.par.best, n = 100,
                 order = 0, seed = NULL, antithetic = TRUE, keep = FALSE,
                 phi = NULL, ...) {
    if (is.numeric(seed))
      set.seed(seed)
    last.par.old <- last.par
    last.par.best.old <- last.par.best
    on.exit({
      last.par <<- last.par.old
      last.par.best <<- last.par.best.old
    })
    h <- spHess(par0, random = TRUE)
    L <- L.created.by.newton
    updateCholesky(L, h)
    rmvnorm <- function(n) {
      u <- matrix(rnorm(ncol(L) * n), ncol(L), n)
      u <- solve(L, u, system = "Lt")
      u <- solve(L, u, system = "Pt")
      as.matrix(u)
    }
    M.5.log2pi <- -0.5 * log(2 * pi)
    logdmvnorm <- function(u) {
      logdetH.5 <- determinant(L, logarithm = TRUE)$modulus
      nrow(h) * M.5.log2pi + logdetH.5 - 0.5 * colSums(u *
                                                         as.matrix(h %*% u))
    }
    eval.target <- function(u, order = 0) {
      par[random] <- u
      f(par, order = order)
    }
    samples <- rmvnorm(n)
    if (antithetic)
      samples <- cbind(samples, -samples)
    log.density.propose <- logdmvnorm(samples)
    samples <- samples + par0[random]
    log.density.target <- -apply(samples, 2, eval.target)
    log.density.target[is.nan(log.density.target)] <- -Inf
    I <- log.density.target - log.density.propose
    M <- max(I)
    if (order >= 1) {
      vec <- exp(I - M)
      p <- vec/sum(vec)
      i <- (p > 0)
      p <- p[i]
      I1 <- apply(samples[, i, drop = FALSE], 2, eval.target,
                  order = 1)[-random, , drop = FALSE]
      gr <- as.vector(I1 %*% p)
      if (order == 1)
        return(gr)
    }
    if (!is.null(phi)) {
      phival <- apply(samples, 2, phi)
      if (is.null(dim(phival)))
        phival <- t(phival)
      p <- exp(I - M)
      p <- p/sum(p)
      ans <- phival %*% p
      return(ans)
    }
    value <- -log(mean(exp(I - M))) - M
    ci <- 1.96 * sd(exp(I - M))/sqrt(n)
    attr(value, "confint") <- -log(mean(exp(I - M)) + c(lower = ci,
                                                        upper = -ci)) - M
    if (keep) {
      attr(value, "samples") <- samples
      attr(value, "nlratio") <- -I
    }
    value
  }
  report <- function(par = last.par) {
    f(par, order = 0, type = "double")
    as.list(reportenv)
  }
  simulate <- function(par = last.par, complete = FALSE) {
    f(par, order = 0, type = "double", do_simulate = TRUE)
    sim <- as.list(reportenv)
    if (complete) {
      ans <- data
      ans[names(sim)] <- sim
    }
    else {
      ans <- sim
    }
    ans
  }
  list(par = par[lfixed()], fn = function(x = last.par[lfixed()],
                                          ...) {
    if (tracepar) {
      cat("par:\n")
      print(x)
    }
    if (!validpar(x)) return(NaN)
    if (is.null(random)) {
      ans <- f(x, order = 0)
      if (!ADreport) {
        if (is.finite(ans) && ans < value.best) {
          last.par.best <<- x
          value.best <<- ans
        }
      }
    } else {
      ans <- try({
        if (MCcontrol$doMC) {
          ff(x, order = 0)
          MC(last.par, n = MCcontrol$n, seed = MCcontrol$seed,
             order = 0)
        } else ff(x, order = 0)
      }, silent = silent)
      if (is.character(ans)) ans <- NaN
    }
    ans
  }, gr = function(x = last.par[lfixed()], ...) {
    if (is.null(random)) {
      ans <- f(x, order = 1)
    } else {
      ans <- try({
        if (MCcontrol$doMC) {
          ff(x, order = 0)
          MC(last.par, n = MCcontrol$n, seed = MCcontrol$seed,
             order = 1)
        } else ff(x, order = 1)
      }, silent = silent)
      if (is.character(ans)) ans <- rep(NaN, length(x))
    }
    if (tracemgc) cat("outer mgc: ", max(abs(ans)), "\n")
    ans
  }, he = function(x = last.par[lfixed()], atomic = usingAtomics()) {
    if (is.null(random)) {
      if (!atomic) return(f(x, order = 2))
      if (is.null(ADGrad)) retape_adgrad()
      return(f(x, type = "ADGrad", order = 1))
    } else {
      stop("Hessian not yet implemented for models with random effects.")
    }
  }, hessian = hessian, method = method, retape = retape, env = env,
  report = report, simulate = simulate, ...)
}
