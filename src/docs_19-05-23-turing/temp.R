#' Temporary plots for presentation

#' Histogram plot
colours <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")

par <- "beta_alpha"
i <- 1

if(!is.null(i)) {
  par_name <- paste0(par, "[", i, "]")

  df_compare <- rbind(
    data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]][i, ])),
    data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]][i, ])),
    data.frame(method = "adam", samples = as.numeric(adam$sample[[par]][i, ])),
    data.frame(method = "tmbstan", samples = as.numeric(tmbstan$mcmc$sample[[par]][i, ]))
  )
} else {
  par_name <- paste0(par)

  df_compare <- rbind(
    data.frame(method = "TMB", samples = as.numeric(tmb$fit$sample[[par]])),
    data.frame(method = "aghq", samples = as.numeric(aghq$quad$sample[[par]])),
    data.frame(method = "adam", samples = as.numeric(adam$sample[[par]])),
    data.frame(method = "tmbstan", samples = as.numeric(tmbstan$mcmc$sample[[par]]))
  )
}

df_compare$method <- factor(df_compare$method, levels = c("TMB", "aghq", "adam", "tmbstan"))

mean <- df_compare %>%
  filter(method == "tmbstan") %>%
  summarise(mean = mean(samples)) %>%
  pull(mean) %>%
  round(digits = 3)

sd <- df_compare %>%
  filter(method == "tmbstan") %>%
  summarise(sd = sd(samples)) %>%
  pull(sd) %>%
  round(digits = 3)

rhat <- tryCatch(round(rhats[[par_name]], 3), error = function(e) return("Missing!"))
ess <- tryCatch(round(ess[[par_name]], 3), error = function(e) return("Missing!"))

histogram_plot <- df_compare %>%
  filter(method != "adam") %>%
  ggplot(aes(x = samples, fill = method, col = method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  facet_grid(method~.) +
  labs(y = "Density", fill = "Method") +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none") +
  labs(
    title = par_name,
    subtitle = paste0("tmbstan: Rhat = ", rhat, ", ESS = ", ess),
    x = ""
  )

ggsave("beta-alpha-histogram.png", histogram_plot, h = 4, w = 6.25)

#' ECDF plot
grid <- seq(from = min(df_compare$samples), to = max(df_compare$samples), length.out = 1000)

tmb_ecdf <- stats::ecdf(filter(df_compare, method == "TMB") %>% pull(samples))
tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "TMB")

aghq_ecdf <- stats::ecdf(filter(df_compare, method == "aghq") %>% pull(samples))
aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "aghq")

adam_ecdf <- stats::ecdf(filter(df_compare, method == "adam") %>% pull(samples))
adam_ecdf_df <- data.frame(x = grid, ecdf = adam_ecdf(grid), method = "adam")

tmbstan_ecdf <- stats::ecdf(filter(df_compare, method == "tmbstan") %>% pull(samples))
tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "tmbstan")

tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
adam_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - adam_ecdf_df$ecdf
tmbstan_ecdf_df$ecdf_diff <- 0

ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
ks_aghq <- absmax(aghq_ecdf_df$ecdf_diff)
ks_adam <- absmax(adam_ecdf_df$ecdf_diff)

ecdf_df <- bind_rows(tmb_ecdf_df, aghq_ecdf_df, adam_ecdf_df, tmbstan_ecdf_df)

ecdf_df$method <- factor(ecdf_df$method, levels = c("TMB", "aghq", "adam", "tmbstan"))

ks_labeller <- function(x) paste0(toString(100 * round(abs(x), 2)), "%")

ecdf_plot <- ecdf_df %>%
  filter(method != "adam") %>%
  ggplot(aes(x = x, y = ecdf, col = method)) +
  geom_line() +
  scale_color_manual(values = colours) +
  labs(x = "", y = "ECDF") +
  guides(col = "none") +
  coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
  theme_minimal() +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

ecdf_diff_plot <- ecdf_df %>%
  filter(method != "adam") %>%
  ggplot(aes(x = x, y = ecdf_diff, col = method)) +
  geom_line() +
  scale_y_continuous(labels = scales::label_percent()) +
  geom_abline(intercept = ks_tmb, slope = 0, col = colours[1], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.15 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = colours[1], alpha = 0.8) +
  geom_abline(intercept = ks_aghq, slope = 0, col = colours[2], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.15 * max(ecdf_df$x), y = ks_aghq, label = ks_labeller(ks_aghq), col = colours[2], alpha = 0.8) +
  scale_color_manual(values = colours) +
  labs(x = "", y = "ECDF difference") +
  guides(col = "none") +
  coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
  theme_minimal() +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

ggsave("ecdf-plot.png", ecdf_plot, h = 4, w = 6.25)
ggsave("ecdf-diff-plot.png", ecdf_diff_plot, h = 4, w = 6.25)

#' Overall KS density plot
ks_summary <- filter(ks_summary, Type != "Hyper")

method1 <- "TMB"
method2 <- "aghq"

ks_method1 <- paste0("KS(", method1, ", tmbstan)")
ks_method2 <- paste0("KS(", method2, ", tmbstan)")

xy_length <- min(1, max(ks_summary[[ks_method1]], ks_summary[[ks_method2]]) + 0.05)

x <- y <- seq(0, xy_length, length.out = 20)

gradient_base <- expand.grid(x, y) %>%
  mutate(diff = Var1 - Var2) %>%
  ggplot(aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = diff), alpha = 0.7) +
  scale_fill_gradientn(
    colours = c("#56B4E9", "white", "#009E73"),
    rescaler = ~ scales::rescale_mid(.x, mid = 0)
  )

densityplot <- gradient_base +
  stat_density_2d(data = ks_summary, aes(x = .data[[ks_method1]], y = .data[[ks_method2]], linetype = Type), col = "black") +
  xlim(0, xy_length) +
  ylim(0, xy_length) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = ks_method1, y = ks_method2) +
  theme_minimal() +
  guides(fill = "none", linetype = "none") +
  theme(
    legend.position = "bottom"
  )

ggsave("ks-density-plot.png", densityplot, h = 4, w = 4)
