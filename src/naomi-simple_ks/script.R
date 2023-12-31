#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_ks")
# setwd("src/naomi-simple_ks")

rmarkdown::render("ks.Rmd")

beta_alpha <- histogram_and_ecdf("beta_alpha", return_df = TRUE)
saveRDS(beta_alpha$df, "beta_alpha.rds")

beta_anc_rho <- histogram_and_ecdf("beta_anc_rho", return_df = TRUE)
saveRDS(beta_anc_rho$df, "beta_anc_rho.rds")

#' Quick plot for thesis!

par <- "ui_lambda_x"
i <- 26

colours <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")

par_name <- paste0(par, "[", i, "]")

df_compare <- rbind(
  data.frame(method = "GEB", samples = as.numeric(tmb$fit$sample[[par]][i, ])),
  data.frame(method = "GPCA-AGHQ", samples = as.numeric(aghq$quad$sample[[par]][i, ])),
  data.frame(method = "NUTS", samples = as.numeric(tmbstan$mcmc$sample[[par]][i, ]))
)

rhat <- signif(rstan::Rhat(tmbstan$mcmc$sample[[par]][i, ]), 3)
ess <- signif(rstan::ess_bulk(tmbstan$mcmc$sample[[par]][i, ]), 3)

df_compare$method <- factor(df_compare$method, levels = c("GEB", "GPCA-AGHQ", "NUTS"))

mean <- signif(mean(filter(df_compare, method == "NUTS")$samples), digits = 3)
sd <- signif(sd(filter(df_compare, method == "NUTS")$samples), digits = 3)

histogram_plot <- ggplot(df_compare, aes(x = samples, fill = method, col = method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  facet_grid(method~.) +
  labs(y = "Density", fill = "Method") +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none") +
  labs(x = "")

grid <- seq(from = min(df_compare$samples), to = max(df_compare$samples), length.out = 1000)

tmb_ecdf <- stats::ecdf(filter(df_compare, method == "GEB") %>% pull(samples))
tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "GEB")

aghq_ecdf <- stats::ecdf(filter(df_compare, method == "GPCA-AGHQ") %>% pull(samples))
aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "GPCA-AGHQ")

tmbstan_ecdf <- stats::ecdf(filter(df_compare, method == "NUTS") %>% pull(samples))
tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "NUTS")

# Add ECDF differences
tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
tmbstan_ecdf_df$ecdf_diff <- 0

ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
ks_aghq <- absmax(aghq_ecdf_df$ecdf_diff)

ecdf_df <- bind_rows(tmb_ecdf_df, aghq_ecdf_df, tmbstan_ecdf_df)

ecdf_df$method <- factor(ecdf_df$method, levels = c("GEB", "GPCA-AGHQ", "NUTS"))

ks_labeller <- function(x) toString(signif(abs(x), 2))

ecdf_plot <- ggplot(ecdf_df, aes(x = x, y = ecdf_diff, col = method)) +
  geom_line() +
  geom_abline(intercept = ks_tmb, slope = 0, col = colours[1], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.4 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = colours[1], alpha = 0.8) +
  geom_abline(intercept = ks_aghq, slope = 0, col = colours[2], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.4 * max(ecdf_df$x), y = ks_aghq, label = ks_labeller(ks_aghq), col = colours[2], alpha = 0.8) +
  scale_color_manual(values = colours) +
  labs(x = "", y = "ECDF difference") +
  guides(col = "none") +
  coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
  theme_minimal() +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

histogram_plot + ecdf_plot

saveRDS(list(par = par_name, rhat = rhat, ess = ess), file = "ui-lambda-x.rds")
ggsave("ui-lambda-x.png", h = 4, w = 6.25)
