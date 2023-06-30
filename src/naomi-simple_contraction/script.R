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
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  annotate("text", x = 36, y = -0.6, size = 3, label = "Prior tighter") +
  annotate("text", x = 36, y = 0.4, size = 3, label = "Posterior tighter") +
  annotate("segment", x = 35, xend = 35, y = 0, yend = -1.0, size = 0.25, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 35, xend = 35, y = 0, yend = 1.0, size = 0.25, arrow = arrow(length = unit(0.2, "cm"))) +
  theme_minimal()

ggsave("posterior-contraction.png", posterior_contraction_plot, h = 6.5, w = 6.25, bg = "white")

#' Which have lower amounts of posterior contraction?
names(subset(posterior_contraction, posterior_contraction < 0.5))
