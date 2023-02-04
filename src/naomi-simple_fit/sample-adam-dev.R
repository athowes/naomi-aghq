#' I've got the fit_adam working with AGHQ k = 1, and know how to sample from the
#' Laplace marginals by using spline interpolation then the inverse CDF sampling
#' method. The next step is simply to bring this together and package it up into
#' a function sample_adam() which outputs in the same format as sample_aghq() say.

compute_pdf_and_cdf <- function(nodes, lps, method = "auto") {
  k <- length(nodes)
  if(k >= 4) method <- "spline"
  if(k < 4) method <- "polynomial"

  rn <- range(nodes)
  rnl <- diff(rn)
  min <- min(rn) - rnl / 2
  max <- max(rn) + rnl / 2

  if(method == "spline") {
    ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
    interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  }

  if(method == "polynomial") {
    interpolant <- as.function(polynom::poly.calc(x = nodes, y = lps))
  }

  finegrid <- seq(from = min, to = max, length.out = 1000)

  df <- data.frame(
    x = finegrid,
    pdf = exp(interpolant(finegrid)),
    cdf = cumsum(exp(interpolant(finegrid))) * c(0, diff(finegrid))
  )

  return(df)
}

sample_cdf <- function(df, M) {
  q <- stats::runif(M)
  samples <- numeric(M)
  for(i in 1:M) samples[i] <- df$x[max(which(df$cdf < q[i]))]
  return(samples)
}

sample_adam <- function(adam, M, verbose = TRUE) {
  if (verbose) print("Sampling from adam")
  r <- adam$quad$obj$env$random
  d <- length(random)
  x_names <- names(adam$quad$obj$env$par[random])
  samps <- matrix(data = NA, nrow = d, ncol = M)
  rownames(samps) <- x_names

  # Mirroring the aghq data structures as in aghq::sample_marginal.marginallaplace
  # Might want to change these as well as in the aghq package eventually

  # Laplace latent field marginals
  for(i in 1:10) {
    marginal <- adam$laplace_marginals[adam$laplace_marginals$index == i, ]
    pdf_and_cdf <- compute_pdf_and_cdf(marginal$x, marginal$lp)
    samps[i, ] <- sample_cdf(pdf_and_cdf, M = M)
  }

  # Hyperparameter marginals
  thetasamples <- list()
  for(j in 1:length(adam$quad$marginals)) {
    marginal <- adam$quad$marginals[[j]]
    colnames(marginal)[grep("theta", colnames(marginal))] <- "theta"
    pdf_and_cdf <- compute_pdf_and_cdf(marginal$theta, marginal$logmargpost)
    thetasamples[[j]] <- unname(sample_cdf(pdf_and_cdf, M = M))
  }

  samp$samps <- samps
  samp$thetasamples <- thetasamples

  # This part replaces samples from TMB with samples from aghq
  if (verbose) print("Rearranging samples")
  smp <- matrix(0, M, length(adam$quad$obj$env$par))
  smp[, r] <- unname(t(samp$samps))
  names(samp$thetasamples) <- names(samp$theta)
  smp[, -r] <- unname(as.matrix(bind_rows(samp$thetasamples)))
  smp <- as.data.frame(smp)
  colnames(smp)[r] <- rownames(samp$samps)
  colnames(smp)[-r] <- names(samp$thetasamples)

  # This part is the same as TMB
  if (verbose) print("Simulating from model")
  sim <- apply(smp, 1, adam$quad$obj$report)
  r <- adam$quad$obj$report()
  if (verbose) print("Returning sample")
  adam$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
  is_vector <- vapply(adam$sample, inherits, logical(1), "numeric")
  adam$sample[is_vector] <- lapply(adam$sample[is_vector], matrix, nrow = 1)
  names(adam$sample) <- names(r)

  adam
}

sample_adam(adam, M = 10)
