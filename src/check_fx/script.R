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

