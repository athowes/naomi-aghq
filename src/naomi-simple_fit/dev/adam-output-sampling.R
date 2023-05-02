#' Debug adam sampling of Naomi outputs
base_grid <- mvQuad::createNIGrid(dim = n_hyper, type = "GHe", level = 1, ndConstruction = "product")
adam <- fit_adam(tmb_inputs_simple, k = 1, basegrid = base_grid)

nsample <- 100

#' Expose sample_adam(adam, M = nsample)
M <- nsample
verbose <- TRUE
if (verbose) print("Sampling from adam")
r <- adam$quad$obj$env$random
d <- length(r)
x_names <- names(adam$quad$obj$env$par[r])
theta_names <- names(adam$quad$obj$env$par[-r])
samps <- matrix(data = NA, nrow = d, ncol = M)
rownames(samps) <- x_names

# Mirroring the aghq data structures as in aghq::sample_marginal.marginallaplace
# Might want to change these as well as in the aghq package eventually

# Laplace latent field marginals
for(i in 1:d) {
  marginal <- adam$laplace_marginals[adam$laplace_marginals$index == i, ]
  # Extra normalisation check in here to avoid any issues with slightly not integrating to one
  pdf_and_cdf <- compute_pdf_and_cdf(marginal$x, marginal$lp_normalised, normalise = TRUE)
  samps[i, ] <- sample_cdf(pdf_and_cdf, M = M)
}

# Hyperparameter marginals
thetasamples <- list()
for(j in 1:length(adam$quad$marginals)) {
  marginal <- adam$quad$marginals[[j]]
  colnames(marginal)[grep("theta", colnames(marginal))] <- "theta"
  # Don't normalise because maybe just one point (for now) and the trapezoid rule will break
  pdf_and_cdf <- compute_pdf_and_cdf(marginal$theta, marginal$logmargpost)
  thetasamples[[j]] <- unname(sample_cdf(pdf_and_cdf, M = M))
}

samp <- list()
samp$samps <- samps
samp$thetasamples <- thetasamples

# This part replaces samples from TMB with samples from aghq
if (verbose) print("Rearranging samples")
smp <- matrix(0, M, length(adam$quad$obj$env$par))
smp[, r] <- unname(t(samp$samps))
names(samp$thetasamples) <- theta_names
smp[, -r] <- unname(as.matrix(bind_rows(samp$thetasamples)))
smp <- as.data.frame(smp)
colnames(smp)[r] <- x_names
colnames(smp)[-r] <- theta_names

# This part is the same as TMB
if (verbose) print("Simulating from model")
sim <- apply(smp, 1, adam$quad$obj$report)
r <- adam$quad$obj$report()
if (verbose) print("Returning sample")
adam$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
is_vector <- vapply(adam$sample, inherits, logical(1), "numeric")
adam$sample[is_vector] <- lapply(adam$sample[is_vector], matrix, nrow = 1)
names(adam$sample) <- names(r)

#' Reproduce histogram from issue showing that the Naomi outputs here are not as they should be
hist(adam$sample$alpha_t1_out[1, ])



