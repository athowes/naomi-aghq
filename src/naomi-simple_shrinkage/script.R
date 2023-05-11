#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_shrinkage")
# setwd("src/naomi-simple_shrinkage")

out <- readRDS("depends/out.rds")
mcmc <- out$mcmc$stanfit

#' Posterior shrinkage is
#' 1 - posterior_variance / prior_variance

#' It's easy to obtain the posterior variance
stan_summary <- rstan::summary(mcmc)$summary
posterior_sd <- stan_summary[, "sd"][1:491] #' Remove the lp__ at the end

#' Obtaining the prior variance is more challenging, but we can do it by modifying the C++ template
TMB::compile("naomi_simple_prior.cpp")
dyn.load(TMB::dynlib("naomi_simple_prior"))
tmbstan::tmbstan()

data_keep <- c("Q_x", "Q_x_rankdef")

obj_prior <- TMB::MakeADFun(
  data = out$inputs$data[data_keep],
  parameters = out$inputs$par_init,
  DLL = "naomi_simple_prior",
  silent = TRUE,
)

prior_mcmc <- tmbstan::tmbstan(obj_prior, chains = 4, iter = 1000)

prior_stan_summary <- rstan::summary(prior_mcmc)$summary
prior_sd <- prior_stan_summary[, "sd"][1:491] #' Remove the lp__ at the end

posterior_shrinkage <- 1 - posterior_sd^2 / prior_sd^2

pdf("posterior-shrinkage.pdf", h = 10, w = 6.25)

data.frame(
  posterior_shrinkage = posterior_shrinkage
) %>%
  tibble::rownames_to_column("par") %>%
  filter(posterior_shrinkage > 0) %>%
  ggplot(aes(y = reorder(par, posterior_shrinkage), x = posterior_shrinkage)) +
    geom_point(alpha = 0.5) +
    labs(x = "Posterior shrinkage", y = "") +
    theme_minimal()

dev.off()

#' Which
names(subset(posterior_shrinkage, posterior_shrinkage < 0.5))
