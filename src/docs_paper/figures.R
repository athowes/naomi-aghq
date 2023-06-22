cols <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Fig A

TMB::compile("2d.cpp")
dyn.load(TMB::dynlib("2d"))

obj <- TMB::MakeADFun(data = list(), parameters = list(theta1 = 0, theta2 = 0), DLL = "2d")

grid <- expand.grid(
  theta1 = seq(-8, 8, length.out = 400),
  theta2 = seq(-8, 8, length.out = 400)
)

ground_truth <- cbind(grid, pdf = apply(grid, 1, function(x) exp(-1 * obj$fn(x))))

opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = 1000, trace = 0)
)

sd_out <- TMB::sdreport(
  obj,
  par.fixed = opt$par,
  getJointPrecision = TRUE
)

mu <- opt$par
cov <- sd_out$cov.fixed

figA0 <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = "lightgrey") +
  coord_fixed(xlim = c(-8, 8), ylim = c(-8, 8), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal() +
  guides(size = "none") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

gg <- mvQuad::createNIGrid(2, "GHe", 3)

add_points <- function(figA0, gg) {

  points <- mvQuad::getNodes(gg) %>%
    as.data.frame() %>%
    mutate(weights = mvQuad::getWeights(gg))

  colnames(points) <- c("theta1", "theta2", "weights")

  figA0 +
    geom_point(
      data = points,
      aes(x = theta1, y = theta2, size = weights),
      alpha = 0.8,
      col = "#009E73",
      inherit.aes = FALSE
    ) +
    scale_size_continuous(range = c(1, 2))
}

figA1 <- add_points(figA0, gg) +
  labs(size = "")

#' Adapt by the mean
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = diag(c(1, 1)), dec.type = 1)

figA2 <- add_points(figA0, gg2) +
  labs(size = "")

#' Adapt by the lower Cholesky
gg3 <- gg
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 2)

figA3 <- add_points(figA0, gg3) +
  labs(size = "")

#' Adapt by the spectral
gg4 <- gg
mvQuad::rescale(gg4, m = mu, C = cov, dec.type = 1)

figA4 <- add_points(figA0, gg4) +
  labs(size = "")

#' PCA-AGHQ
gg5 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg5, m = mu, C = cov, dec.type = 1)

figA5 <- add_points(figA0, gg5) +
  labs(size = "")

figA <- (figA1 + figA2) / (figA3 + figA4 + figA5) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position  = c(0.15, 0.95))

#' Fig B

df_compare <- readRDS("depends/beta_alpha.rds") %>%
  filter(method != "adam")

mean <- df_compare %>%
  filter(method == "tmbstan") %>%
  summarise(mean = mean(samples)) %>%
  pull(mean) %>%
  signif(digits = 2)

sd <- df_compare %>%
  filter(method == "tmbstan") %>%
  summarise(sd = sd(samples)) %>%
  pull(sd) %>%
  signif(digits = 2)

histogram <- df_compare %>%
  mutate(method = fct_recode(method, "PCA-AGHQ" = "aghq", "NUTS"= "tmbstan")) %>%
  ggplot(aes(x = samples, fill = method, col = method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  facet_grid(method~.) +
  labs(x = "beta_alpha", y = "Density", fill = "Method") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none")

grid <- seq(from = min(df_compare$samples), to = max(df_compare$samples), length.out = 1000)

tmb_ecdf <- stats::ecdf(filter(df_compare, method == "TMB") %>% pull(samples))
tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "TMB")

aghq_ecdf <- stats::ecdf(filter(df_compare, method == "aghq") %>% pull(samples))
aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "aghq")

tmbstan_ecdf <- stats::ecdf(filter(df_compare, method == "tmbstan") %>% pull(samples))
tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "tmbstan")

# Add ECDF differences
tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
tmbstan_ecdf_df$ecdf_diff <- 0

absmax <- function(x) x[which.max(abs(x))]

ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
ks_aghq <- absmax(aghq_ecdf_df$ecdf_diff)

ecdf_df <- bind_rows(tmb_ecdf_df, aghq_ecdf_df, tmbstan_ecdf_df)

ecdf_df$method <- factor(ecdf_df$method, levels = c("TMB", "aghq", "tmbstan"))

ks_labeller <- function(x) toString(signif(abs(x), 2))

ecdf_diff <- ggplot(ecdf_df, aes(x = x, y = ecdf_diff, col = method)) +
  geom_line() +
  geom_abline(intercept = ks_tmb, slope = 0, col = cols[1], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = cols[1], alpha = 0.8) +
  geom_abline(intercept = ks_aghq, slope = 0, col = cols[2], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_aghq, label = ks_labeller(ks_aghq), col = cols[2], alpha = 0.8) +
  scale_color_manual(values = cols) +
  labs(x = "beta_alpha", y = "ECDF difference") +
  guides(col = "none") +
  coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
  theme_minimal() +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

figB <- histogram + ecdf_diff
