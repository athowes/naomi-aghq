#' Debugging tmbstan not producing output samples properly

#' Start of if(tmbstan) loop

#' A small number of iterations
niter <- 100

mcmc <- fit_tmbstan(tmb_inputs, chains = 4, iter = niter, thin = nthin, cores = 4)

#' Expose sample_tmbstan(mcmc, M = NULL, verbose = TRUE)

mcmc <- mcmc
M <- NULL
verbose <- TRUE

if (verbose) print("Samples already available for tmbstan.")
if (verbose) print("M ignored. Instead length of MCMC chains used.")
if (verbose) print("Simulating from model")
smp <- as.matrix(mcmc$stanfit)
smp <- smp[, colnames(smp) != "lp__"] # Remove the lp__ column
colnames(smp) <- names(mcmc$obj$env$par)
sim <- apply(smp, 1, mcmc$obj$report)
r <- mcmc$obj$report()
if (verbose) print("Returning sample")
sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
is_vector <- vapply(sample, inherits, logical(1), "numeric")
sample[is_vector] <- lapply(sample[is_vector], matrix, nrow = 1)
names(sample) <- names(r)
mcmc$sample <- sample
mcmc
