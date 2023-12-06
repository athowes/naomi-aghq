#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_contraction")
# setwd("src/naomi-simple_contraction")

tmb <- readRDS("depends/tmb.rds")
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

#' Takes a while
prior_mcmc <- tmbstan::tmbstan(obj_prior, chains = 4, iter = 20000)

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
  labs(x = "", y = "Posterior contraction", col = "", shape = "") +
  scale_color_manual(values = c("#56B4E9", "#009E73")) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25, col = "grey40") +
  annotate("text", x = 46, y = -0.8, size = 3, label = "Prior tighter", col = "grey40") +
  annotate("text", x = 46, y = 0.4, size = 3, label = "Posterior tighter", col = "grey40") +
  theme_minimal()

ggsave("posterior-contraction.png", posterior_contraction_plot, h = 6.5, w = 6.25, bg = "white")

#' Which have lower amounts of posterior contraction?
names(subset(posterior_contraction, posterior_contraction < 0.5))

#' Calculate the prior standard deviations manually for the hyperparameters
nsim <- 100000
phi_A <- rbeta(nsim, 0.5, 0.5)
logit_phi_A_sd <- sd(qlogis(phi_A))

logit_phi_B <- rnorm(nsim, 0, 2.582)
logit_phi_B_sd <- sd(logit_phi_B)

sigma_A <- abs(rnorm(nsim, 0, 2.5))
log_sigma_A_sd <- sd(log(sigma_A))

sigma_B <- abs(rnorm(nsim, 0, 1.0))
log_sigma_B_sd <- sd(log(sigma_B))

manual_hyper_prior_sd <- list(
  "logit_phi_rho_x" = logit_phi_A_sd,
  "log_sigma_rho_x" = log_sigma_A_sd,
  "logit_phi_rho_xs" = logit_phi_A_sd,
  "log_sigma_rho_xs" = log_sigma_A_sd,
  "logit_phi_rho_a" = logit_phi_B_sd,
  "log_sigma_rho_a" = log_sigma_A_sd,
  "logit_phi_rho_as" = logit_phi_B_sd,
  "log_sigma_rho_as" = log_sigma_A_sd,
  "log_sigma_rho_xa" = log_sigma_A_sd,
  "logit_phi_alpha_x" = logit_phi_A_sd,
  "log_sigma_alpha_x" = log_sigma_A_sd,
  "logit_phi_alpha_xs" = logit_phi_A_sd,
  "log_sigma_alpha_xs" = log_sigma_A_sd,
  "logit_phi_alpha_a" = logit_phi_B_sd,
  "log_sigma_alpha_a" = log_sigma_A_sd,
  "logit_phi_alpha_as" = logit_phi_B_sd,
  "log_sigma_alpha_as" = log_sigma_A_sd,
  "log_sigma_alpha_xa" = log_sigma_A_sd,
  "OmegaT_raw" = 1.0,
  "log_betaT" = log_sigma_B_sd,
  "log_sigma_lambda_x" = log_sigma_A_sd,
  "log_sigma_ancrho_x" = log_sigma_A_sd,
  "log_sigma_ancalpha_x" = log_sigma_A_sd,
  "log_sigma_or_gamma" = log_sigma_A_sd
)

nuts_manual_comparison_df <- prior_sd[names(tmb$fit$obj$par)] %>%
  data.frame() %>%
  tibble::rownames_to_column("par") %>%
  rename("NUTS" = ".") %>%
  mutate(Direct = unlist(manual_hyper_prior_sd))

nuts_manual_comparison_df %>%
  pivot_longer(
    cols = c("NUTS", "Direct"),
    names_to = "method",
    values_to = "sd"
  ) %>%
  ggplot(aes(x = par, y = sd, col = method, shape = method)) +
    geom_point(size = 2) +
    coord_flip() +
    scale_color_manual(values = c("#E69F00", "#F0E442")) +
    labs(x = "Hyperparameter", y = "Prior standard deviation", col = "Calculation\nmethod", shape = "Calculation\nmethod") +
    theme_minimal()

ggsave("nuts-hand-comparison.png", h = 4.5, w = 6.25, bg = "white")
