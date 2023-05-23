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
  mutate(type = "Hyper")

df_latent <- data.frame(
  posterior_contraction = posterior_contraction
) %>%
  tibble::rownames_to_column("par") %>%
  filter(!(par %in% names(tmb$fit$obj$par))) %>%
  mutate(par_name = ifelse(str_detect(par, "\\["), str_extract(par, ".*(?=\\[)"), par)) %>%
  group_by(par_name) %>%
  summarise(posterior_contraction = mean(posterior_contraction)) %>%
  mutate(type = "Latent (averaged)") %>%
  rename(par = par_name)

pdf("posterior-contraction.pdf", h = 5, w = 8)

bind_rows(df_hyper, df_latent) %>%
  ggplot(aes(y = reorder(par, posterior_contraction), x = posterior_contraction)) +
  facet_grid(cols = vars(type), scales = "free_x", ) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  labs(x = "Posterior contraction", y = "") +
  scale_x_continuous(limits = c(-0.3, 1), breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1))

dev.off()

#' Which have lower amounts of posterior contraction?
names(subset(posterior_contraction, posterior_contraction < 0.5))

#' How does the size of standard deviations compare in TMB to AGHQ?
tmb_sd <- lapply(tmb$fit$sample[df_latent$par], function(row) data.frame("TMB" = matrixStats::rowSds(row))) %>%
  dplyr::bind_rows(.id = "par")

aghq_sd <- lapply(aghq$quad$sample[df_latent$par], function(row) data.frame("aghq" = matrixStats::rowSds(row))) %>%
  dplyr::bind_rows(.id = "par")

tmbstan_sd <- lapply(tmbstan$mcmc$sample[df_latent$par], function(row) data.frame("tmbstan" = matrixStats::rowSds(row))) %>%
  dplyr::bind_rows(.id = "par")

df_sd <- tmb_sd
df_sd$aghq <- aghq_sd$aghq
df_sd$tmbstan <- tmbstan_sd$tmbstan

pdf("sd-comparison.pdf", h = 5, w = 6.25)

plot_aghq_sd <- ggplot(df_sd, aes(x = tmbstan, y = aghq)) +
  geom_point(alpha = 0.5) +
  labs(x = "SD (tmbstan)", y = "SD (aghq)") +
  geom_abline(slope = 1, intercept = 0, col = "#CC79A7", linetype = "dashed") +
  theme_minimal()

plot_tmb_sd <- ggplot(df_sd, aes(x = tmbstan, y = TMB)) +
  geom_point(alpha = 0.5) +
  labs(x = "SD (tmbstan)", y = "SD (TMB)") +
  geom_abline(slope = 1, intercept = 0, col = "#CC79A7", linetype = "dashed") +
  theme_minimal()

plot_tmb_sd + plot_aghq_sd

dev.off()
