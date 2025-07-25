library(tidyverse)
library(patchwork)

cols <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Fig A

TMB::compile("2d.cpp")
dyn.load(TMB::dynlib("2d"))

obj <- TMB::MakeADFun(data = list(), parameters = list(theta1 = 0, theta2 = 0), DLL = "2d")

box_lower <- -5
box_upper <- 10
box_size <- box_upper - box_lower

grid <- expand.grid(
  theta1 = seq(box_lower, box_upper, length.out = box_size * 50),
  theta2 = seq(box_lower, box_upper, length.out = box_size * 50)
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
  coord_fixed(xlim = c(box_lower, box_upper), ylim = c(box_lower, box_upper), ratio = 1) +
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
    scale_size_continuous(range = c(1, 2)) +
    theme(
      plot.caption = element_text(hjust = 0)
    )
}

figA1 <- add_points(figA0, gg) +
  labs(subtitle = "A", size = "", caption = "GHQ")

#' Adapt by the spectral
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = cov, dec.type = 1)

figA2 <- add_points(figA0, gg2) +
  labs(subtitle = "B", size = "", caption = "AGHQ")

#' PCA-AGHQ
gg3 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 1)

lambda <- eigen(cov)$values
cumsum(lambda) / sum(lambda)

xstart <- 6.2
ystart <- -2.3

x1end <- xstart + 4 * eigen(cov)$vectors[1, 1]
y1end <- ystart + 4 * eigen(cov)$vectors[2, 1]

x2end <- xstart + 1 * eigen(cov)$vectors[1, 2]
y2end <- ystart + 1 * eigen(cov)$vectors[2, 2]

figA3 <- add_points(figA0, gg3) +
  geom_segment(aes(x = xstart, y = ystart, xend = x1end, yend = y1end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x1end + 1, y = y1end - 3, label = "95%", col = "darkgrey") +
  geom_segment(aes(x = xstart, y = ystart, xend = x2end, yend = y2end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x2end, y = y2end - 2, label = "5%", col = "darkgrey") +
  labs(subtitle = "C", size = "", caption = "PCA-AGHQ")

cowplot::plot_grid(figA1, figA2, figA3, ncol = 3)

ggsave("figA.png", h = 3, w = 6.25)

#' Fig B

# aghq <- readRDS("depends/aghq.rds")
tmb <- readRDS("depends/tmb.rds")
outputs <- tmb$outputs

#' Adapted from figure from Naomi paper (Eaton 2021)
#' https://github.com/jeffeaton/naomi-model-paper/blob/master/analysis.R

indicators <- naomi::add_output_labels(outputs) %>%
  left_join(
    select(outputs$meta_area, area_level, area_id, center_x, center_y),
    by = c("area_level", "area_id")
  ) %>%
  sf::st_as_sf() %>%
  mutate(
    area_level_label = fct_reorder(area_level_label, area_level),
    age_group_label = fct_reorder(age_group_label, as.integer(factor(age_group)))
  )

figB1data <- indicators %>%
  filter(
    indicator == "prevalence",
    age_group == "Y015_049",
    sex == "both",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  )

figB1 <- figB1data %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  scale_fill_viridis_c(
    option = "C", direction = -1,
    begin = 0.1, end = 0.9,
    labels = scales::label_percent(1)
  ) +
  expand_limits(fill = 0) +
  labs(subtitle = "A", fill = "", caption = "HIV prevalence") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8), angle = 90),
    legend.text = element_text(size = rel(0.8)),
    # legend.position = "bottom",
    legend.direction = "vertical",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(1, "line")
  )

figB2data <- indicators %>%
  filter(
    indicator == "art_coverage",
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    sex == "both",
    area_level == 4
  )

figB2 <- figB2data %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  scale_fill_viridis_c(
    option = "D", direction = -1, begin = 0.05, end = 0.9,
    labels = scales::label_percent(1), limits = c(0.6, 0.853)
  ) +
  labs(subtitle = "B", fill = "", caption = "ART coverage") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8), angle = 90),
    legend.text = element_text(size = rel(0.8)),
    # legend.position = "bottom",
    legend.direction = "vertical",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(1, "line")
  )

figB3data <- indicators %>%
  filter(
    indicator %in% c("infections", "incidence"),
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    sex == "both",
    area_level == 4
  )

figB3 <- figB3data %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_wider(c(area_id, area_name), names_from = indicator, values_from = mean) %>%
  left_join(
    select(outputs$meta_area, area_id, center_x, center_y),
    by = "area_id"
  ) %>%
  sf::st_as_sf() %>%
  ggplot(aes(x = center_x, y = center_y, color = incidence, size = infections)) +
  geom_sf(size = 0.1, color = "grey30") +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(
    option = "B", direction = -1,
    begin = 0.05, end = 0.9,
    labels = scales::label_number(scale = 1000)
  ) +
  scale_size_area(max_size = 8) +
  labs(
    subtitle = "C",
    caption = "New HIV cases and HIV incidence",
    col = "",
    size = "",
    x = element_blank(),
    y = element_blank()
  ) +
  expand_limits(color = 0) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8), angle = 90),
    legend.text = element_text(size = rel(0.8)),
    # legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "vertical",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(1, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-10, -10, -10, -10)
  )

ggsave(
  "figB.png",
  plot = figB1 + figB2 + figB3,
  width = 6.25,
  height = 3
)

#' Fig C

df_compare <- readRDS("depends/beta_alpha.rds")

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

histogram + ecdf_diff

ggsave("figC.png", h = 4, w = 6.25, background = "white")

#' Abstract figure
top_row <- cowplot::plot_grid(figA1 + labs(subtitle = ""), figA2 + labs(subtitle = ""), figA3 + labs(subtitle = ""), ncol = 3)
bottom_row <- {figB1 + labs(subtitle = "") + figB2 + labs(subtitle = "") + figB3 + labs(subtitle = "")}
ggsave("abstract.png", cowplot::plot_grid(top_row, bottom_row, nrow = 2), h = 6, w = 6.25, bg = "white")
