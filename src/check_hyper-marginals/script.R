#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_hyper-marginals")
# setwd("src/check_hyper-marginals")

tmbstan <- readRDS("depends/tmbstan.rds")

hyper <- names(tmbstan$mcmc$obj$par)
stopifnot(length(hyper) == 24)

pdf("hist.pdf", h = 10, w = 10)

bayesplot::mcmc_hist(tmbstan$mcmc$stanfit, pars = hyper) +
  theme_minimal()

dev.off()

pdf("pairs.pdf", h = 10, w = 10)

bayesplot::mcmc_pairs(tmbstan$mcmc$stanfit, pars = hyper[1:6], off_diag_fun = "hex")
bayesplot::mcmc_pairs(tmbstan$mcmc$stanfit, pars = hyper[7:12], off_diag_fun = "hex")
bayesplot::mcmc_pairs(tmbstan$mcmc$stanfit, pars = hyper[13:18], off_diag_fun = "hex")
bayesplot::mcmc_pairs(tmbstan$mcmc$stanfit, pars = hyper[19:24], off_diag_fun = "hex")

dev.off()
