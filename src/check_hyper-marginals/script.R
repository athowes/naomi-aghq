#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_hyper-marginals")
# setwd("src/check_hyper-marginals")

tmbstan <- readRDS("depends/tmbstan.rds")
aghq <- readRDS("depends/aghq.rds")

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

#' So it looks like things are not totally Gaussian -- that's good
#' Now look at the positions of the aghq nodes on here

nodes_samples_comparison <- function(par) {
  tmbstan_samples <- as.numeric(unlist(tmbstan$mcmc$sample[par]))
  aghq_nodes <- aghq$quad$modesandhessians[[par]]

  l <- floor(min(tmbstan_samples))
  u <- ceiling(max(tmbstan_samples))

  plot_nodes <- data.frame(index = 1:length(aghq_nodes), node = aghq_nodes) %>%
    ggplot(aes(y = index, x = node)) +
    geom_point(alpha = 0.2) +
    xlim(l, u) +
    labs(y = "Index of node", x = paste0(par), subtitle = "AGHQ nodes") +
    theme_minimal()

  plot_samples <- data.frame(x = tmbstan_samples) %>%
    ggplot(aes(x = x)) +
    geom_histogram(fill = "lightgrey", col = "darkgrey") +
    xlim(l, u) +
    labs(x = paste0(par), y = "", subtitle = "HMC samples") +
    theme_minimal()

  plot_nodes / plot_samples
}

pdf("nodes-samples-comparison.pdf", h = 7, w = 6.25)

lapply(1:length(hyper), function(i) nodes_samples_comparison(hyper[i]))

dev.off()
