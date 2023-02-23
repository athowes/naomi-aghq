#' Want to calculate the ESS of indicators produced from sample_method functions

out <- readRDS("depends/out.rds")
mcmc <- out$mcmc

#' Check to make sure that this MCMC has the samples in it
stopifnot(!is.null(mcmc$sample))

#' Which parameters are there in `mcmc$sample` that are not sampled with MCMC?
names(mcmc$sample)[!(names(mcmc$sample) %in% unique(names(mcmc$obj$env$par)))]

#' Pick a parameter in the MCMC
x <- mcmc$sample$log_sigma_rho_x
str(x)
length(x)
nrow(x)
ncol(x)

#' Test that ess_****() gets the same thing as in the rstan summary
rstan::ess_bulk(t(x))
rstan::ess_tail(t(x))

rstan::summary(mcmc$stanfit)$summary["log_sigma_rho_x", ]
rstan::summary(mcmc$stanfit)$c_summary

bayesplot::mcmc_trace(mcmc$stanfit, pars = "log_sigma_rho_x")

