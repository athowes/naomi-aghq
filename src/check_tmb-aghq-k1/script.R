#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_tmb-aghq-k1")
# setwd("src/check_tmb-aghq-k1")

tmb <- readRDS("depends/tmb.rds")
aghqk1 <- readRDS("depends/aghqk1.rds")
aghq <- readRDS("depends/aghq.rds")

par_samp_matrix <- function(sample) {
  x <- sample[!(stringr::str_ends(names(sample), "_out") | stringr::str_ends(names(sample), "_ll"))]
  do.call(rbind, x)
}

tmb_samples <- par_samp_matrix(tmb$fit$sample)
aghqk1_samples <- par_samp_matrix(aghqk1$quad$sample)
aghq_samples <- par_samp_matrix(aghq$quad$sample)

n <- nrow(tmb_samples)
ks_results <- lapply(1:n, function(i) ks.test(tmb_samples[i, ], aghqk1_samples[i, ]))
p_values <- sapply(ks_results, function(x) x$p.value)

#' Expect to show here that TMB and AGHQ k = 1 are not significantly different
data.frame(p_values = p_values) %>%
  filter(p_values != 1) %>%
  ggplot(aes(x = p_values)) +
    geom_histogram(fill = "lightgrey", col = "grey") +
    labs(x = "p-value", y = "Count", title = "KS test p-values for TMB against AGHQ with k = 1") +
    theme_minimal()

#' Expect to show here that TMB and AGHQ k = 3, s = 8 are significantly different
ks_results_aghq <- lapply(1:n, function(i) ks.test(tmb_samples[i, ], aghq_samples[i, ]))
p_values_aghq  <- sapply(ks_results_aghq , function(x) x$p.value)

data.frame(p_values = p_values_aghq) %>%
  filter(p_values != 1) %>%
  ggplot(aes(x = p_values)) +
  geom_histogram(fill = "lightgrey", col = "grey") +
  labs(x = "p-value", y = "Count", title = "KS test p-values for TMB against AGHQ with k = 3, s = 8") +
  theme_minimal()
