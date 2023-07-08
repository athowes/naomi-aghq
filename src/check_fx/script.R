#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_fx")
# setwd("src/check_fx")

tmb <- readRDS("depends/tmb.rds")
aghq <- readRDS("depends/aghq.rds")
tmbstan <- readRDS("depends/tmbstan.rds")

latent_pars <- unique(names(tmb$fit$par.full))[!(unique(names(tmb$fit$par.full)) %in% names(tmb$fit$par))]
output_pars <- c("rho_t1_out", "lambda_t1_out", "alpha_t1_out")

#' Start with f where we expect performance in x to translate into performance in f(x)
#' Chose f(x) = sum(x)

latent_sum <- function(samples) {
  x <- samples[latent_pars]
  x <- do.call(rbind, x)
  colSums(x)
}

tmb_sum <- latent_sum(tmb$fit$sample)
aghq_sum <- latent_sum(aghq$quad$sample)
tmbstan_sum <- latent_sum(tmbstan$mcmc$sample)

df <- bind_rows(
  data.frame(x = tmb_sum, method = "TMB"),
  data.frame(x = aghq_sum, method = "PCA-AGHQ"),
  data.frame(x = tail(tmbstan_sum, n = 1000), method = "NUTS")
)

ggplot(df, aes(x = x, y = ..density.., fill = method)) +
  geom_histogram() +
  facet_grid(~ method) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#56B4E9")) +
  labs(x = "Sum of latent field", y = "Density", fill = "") +
  theme_minimal()

ggsave("sum-latent.png", h = 4, w = 6.25, bg = "white")

df %>%
  group_by(method) %>%
  summarise(
    mean = mean(x),
    sd = sd(x)
  )

#' Now reproduce HIV output from scratch. This may be useful:
#' https://github.com/mrc-ide/naomi/blob/master/R/tmb-model-r-implementation.R

calculate_rho_t1 <- function(p, d) {
  mu_rho <- d$X_rho %*% p$beta_rho +
    d$logit_rho_offset +
    d$Z_rho_x %*% p$u_rho_x +
    d$Z_rho_xs %*% p$u_rho_xs +
    d$Z_rho_a %*% p$u_rho_a +
    d$Z_rho_as %*% p$u_rho_as

  mu_rho <- as.vector(mu_rho)

  rho_15to49f_t1 <- d$X_15to49f %*% (plogis(mu_rho) * d$population_t1) / (d$X_15to49f %*% d$population_t1)
  mu_rho_paed <- d$X_paed_rho_ratio %*% rho_15to49f_t1 + d$paed_rho_ratio_offset
  mu_rho_paed <- as.vector(mu_rho_paed)
  mu_rho_paed <- qlogis(mu_rho_paed)
  mu_rho <- mu_rho + mu_rho_paed

  rho_t1 <- plogis(mu_rho)
  return(rho_t1)
}

f <- function(sample) {
  d <- tmb$inputs$data
  n <- ncol(sample$beta_rho)
  X <- matrix(nrow = n, ncol = 1088)
  for(i in 1:n) {
    p <- lapply(sample[latent_pars], function(x) x[, i])
    X[i, ] <- calculate_rho_t1(p, d)
  }
  return(X)
}

X_tmb <- f(tmb$fit$sample)
X_aghq <- f(aghq$quad$sample)
X_tmbstan <- f(tmbstan$mcmc$sample)

tmb_error <- colMeans(X_tmb) - colMeans(X_tmbstan)
aghq_error <- colMeans(X_aghq) - colMeans(X_tmbstan)

sqrt(mean(tmb_error^2))
sqrt(mean(aghq_error^2))

tmb_error <- apply(X_tmb, 2, sd) - apply(X_tmbstan, 2, sd)
aghq_error <- apply(X_aghq, 2, sd) - apply(X_tmbstan, 2, sd)

plot(tmb_error)
plot(aghq_error)

sqrt(mean(tmb_error^2))
sqrt(mean(aghq_error^2))
