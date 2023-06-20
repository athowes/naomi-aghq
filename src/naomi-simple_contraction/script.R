#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_contraction")
# setwd("src/naomi-simple_contraction")

tmb <- readRDS("depends/tmb.rds")
aghq <- readRDS("depends/aghq.rds")
tmbstan <- readRDS("depends/tmbstan.rds")
mcmc <- tmbstan$mcmc$stanfit

#' Posterior contraction is
#' 1 - posterior_variance / prior_variance

#' It's easy to obtain the posterior variance
stan_summary <- rstan::summary(mcmc)$summary
posterior_sd <- stan_summary[, "sd"][1:491] #' Remove the lp__ at the end

#' Obtaining the prior variance is more challenging, but we can do it by modifying the C++ template
TMB::compile("naomi_simple_prior.cpp")
dyn.load(TMB::dynlib("naomi_simple_prior"))

data_keep <- c("Q_x", "Q_x_rankdef")

obj_prior <- TMB::MakeADFun(
  data = tmbstan$inputs$data[data_keep],
  parameters = tmbstan$inputs$par_init,
  DLL = "naomi_simple_prior",
  silent = TRUE,
)

prior_mcmc <- tmbstan::tmbstan(obj_prior, chains = 4, iter = 1000)

prior_stan_summary <- rstan::summary(prior_mcmc)$summary
prior_sd <- prior_stan_summary[, "sd"][1:491] #' Remove the lp__ at the end

posterior_contraction <- 1 - posterior_sd^2 / prior_sd^2

df_hyper <- data.frame(
  posterior_contraction = posterior_contraction
) %>%
  tibble::rownames_to_column("par") %>%
  filter(par %in% names(tmb$fit$obj$par)) %>%
  mutate(type = "Hyperparameter")

df_latent <- data.frame(
  posterior_contraction = posterior_contraction
) %>%
  tibble::rownames_to_column("par") %>%
  filter(!(par %in% names(tmb$fit$obj$par))) %>%
  mutate(par_name = ifelse(str_detect(par, "\\["), str_extract(par, ".*(?=\\[)"), par)) %>%
  group_by(par_name) %>%
  summarise(posterior_contraction = mean(posterior_contraction)) %>%
  mutate(type = "Latent field") %>%
  rename(par = par_name)

posterior_contraction_plot <- bind_rows(df_hyper, df_latent) %>%
  ggplot(aes(x = reorder(par, posterior_contraction), y = posterior_contraction, col = type, shape = type)) +
  geom_point() +
  coord_flip() +
  labs(x = "", y = "Posterior contraction", col = "Type", shape = "Type") +
  scale_color_manual(values = c("#56B4E9", "#009E73")) +
  scale_y_continuous(limits = c(-0.3, 1), breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal()

ggsave("posterior-contraction.png", posterior_contraction_plot, h = 6.5, w = 6.25)

#' Which have lower amounts of posterior contraction?
names(subset(posterior_contraction, posterior_contraction < 0.5))

#' How does the size of standard deviations compare in TMB to AGHQ?
df_tmb <- lapply(tmb$fit$sample[df_latent$par], function(row) {
  data.frame("sd" = matrixStats::rowSds(row), "mean" = rowMeans(row)) %>%
    tibble::rowid_to_column("id")
  }) %>%
  dplyr::bind_rows(.id = "par") %>%
  mutate(method = "TMB")

df_aghq <- lapply(aghq$quad$sample[df_latent$par], function(row) {
  data.frame("sd" = matrixStats::rowSds(row), "mean" = rowMeans(row)) %>%
    tibble::rowid_to_column("id")
  }) %>%
  dplyr::bind_rows(.id = "par") %>%
  mutate(method = "PCA-AGHQ")

df_tmbstan <- lapply(tmbstan$mcmc$sample[df_latent$par], function(row) {
  data.frame("sd" = matrixStats::rowSds(row), "mean" = rowMeans(row)) %>%
    tibble::rowid_to_column("id")
  }) %>%
  dplyr::bind_rows(.id = "par") %>%
  mutate(method = "NUTS")

df <- bind_rows(df_tmb, df_aghq, df_tmbstan)

df_plot <- df %>%
  pivot_longer(cols = c("sd", "mean"), names_to = "indicator", values_to = "estimate") %>%
  pivot_wider(values_from = "estimate", names_from = "method") %>%
  pivot_longer(cols = c("TMB", "PCA-AGHQ"), names_to = "method", values_to = "approximate") %>%
  rename("truth" = "NUTS") %>%
  mutate(
    indicator = fct_recode(indicator, "Posterior mean estimate" = "mean", "Posterior SD estimate" = "sd"),
    method = fct_relevel(method, "TMB", "PCA-AGHQ")
  )

df_metrics <- df_plot %>%
  group_by(method, indicator) %>%
  summarise(
    rmse = sqrt(mean((truth - approximate)^2)),
    mae = mean(abs(truth - approximate))
  )

df_metrics_pct <- df_metrics %>%
  ungroup() %>%
  group_by(indicator) %>%
  summarise(
    rmse_diff = 100 * diff(rmse) / max(rmse),
    mae_diff = 100 * diff(mae) / max(mae)
  )

df_metrics <- df_metrics %>%
  mutate(
    label = ifelse(
      method == "PCA-AGHQ",
      paste0("RMSE: ", round(rmse, 2), " (", round(df_metrics_pct$rmse_diff), "%)", "\nMAE: ", round(mae, 2), " (", round(df_metrics_pct$mae_diff), "%)"),
      paste0("RMSE: ", round(rmse, 2), "\nMAE: ", round(mae, 2))
    )
  )

write_csv(df_metrics, "mean-sd.csv")

mean_sd_plot <- ggplot(df_plot, aes(x = truth, y = approximate)) +
  geom_point(shape = 1, alpha = 0.4) +
  facet_grid(indicator ~ method) +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = df_metrics, aes(x = -Inf, y = Inf, label = label), size = 3, hjust = 0, vjust = 1.5) +
  labs(x = "NUTS", y = "") +
  theme_minimal()

ggsave("mean-sd.png", mean_sd_plot, h = 6, w = 6.25)

#' Split into two plots for presentations etc.

mean_plot <- df_plot %>%
  filter(indicator == "Posterior mean estimate") %>%
  ggplot(aes(x = truth, y = approximate)) +
  geom_point(shape = 1, alpha = 0.4) +
  facet_grid(indicator ~ method) +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(
    data = filter(df_metrics, indicator == "Posterior mean estimate"),
    aes(x = -Inf, y = Inf, label = label), size = 3, hjust = 0, vjust = 1.5
  ) +
  labs(x = "NUTS", y = "") +
  theme_minimal()

ggsave("mean.png", mean_plot, h = 4, w = 6.25)

sd_plot <- df_plot %>%
  filter(indicator == "Posterior SD estimate") %>%
  ggplot(aes(x = truth, y = approximate)) +
  geom_point(shape = 1, alpha = 0.4) +
  facet_grid(indicator ~ method) +
  coord_fixed(ratio = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(
    data = filter(df_metrics, indicator == "Posterior SD estimate"),
    aes(x = -Inf, y = Inf, label = label), size = 3, hjust = 0, vjust = 1.5
  ) +
  labs(x = "NUTS", y = "") +
  theme_minimal()

ggsave("sd.png", sd_plot, h = 4, w = 6.25)
