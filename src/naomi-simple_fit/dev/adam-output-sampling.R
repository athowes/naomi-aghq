#' Debug adam sampling of Naomi outputs
base_grid <- mvQuad::createNIGrid(dim = n_hyper, type = "GHe", level = 1, ndConstruction = "product")
adam <- fit_adam(tmb_inputs_simple, k = 1, basegrid = base_grid)

nsample <- 1000

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

#' Is it that the latent field and hyper samples are wrong, or that the method for sampling outputs is wrong?
quad <- fit_aghq(tmb_inputs_simple, k = 1)

#' Expose sample_aghq(quad, M = nsample)

# Note that with k = 1, sample_marginal just returns the mode for the hypers
# This needs to be debugged, maybe?
if (verbose) print("Sampling from aghq")
samp <- aghq::sample_marginal(quad, M)

# This part replaces samples from TMB with samples from aghq
if (verbose) print("Rearranging samples")
r <- quad$obj$env$random
smp <- matrix(0, M, length(quad$obj$env$par))
smp[, r] <- unname(t(samp$samps))
names(samp$thetasamples) <- names(samp$theta)
smp[, -r] <- unname(as.matrix(bind_rows(samp$thetasamples)))
smp <- as.data.frame(smp)
colnames(smp)[r] <- rownames(samp$samps)
colnames(smp)[-r] <- names(samp$thetasamples)

# This part is the same as TMB
if (verbose) print("Simulating from model")
sim <- apply(smp, 1, quad$obj$report)
r <- quad$obj$report()
if (verbose) print("Returning sample")
quad$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
is_vector <- vapply(quad$sample, inherits, logical(1), "numeric")
quad$sample[is_vector] <- lapply(quad$sample[is_vector], matrix, nrow = 1)
names(quad$sample) <- names(r)

#' This does look different
hist(quad$sample$alpha_t1_out[1, ])

#' The problem is already there at before apply(smp, 1, obj$report)
#' It is also already there before if (verbose) print("Rearranging samples")
#' Must be during the sampling of draws from adam

#' The hyperparameter samples look the same
lapply(samp$thetasamples, head)
lapply(thetasamples, head)

#' So it must be a difference in the latent field parameters?
hist(samps[1, ])
hist(samp$samps[1, ])

hist(samps[2, ])
hist(samp$samps[2, ])
