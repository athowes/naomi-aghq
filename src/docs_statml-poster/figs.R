#' Figure 1

#' A screenshot from https://naomi.unaids.org/

#' Figure 2

TMB::compile("2d.cpp")
dyn.load(TMB::dynlib("2d"))

obj <- TMB::MakeADFun(data = list(), parameters = list(theta1 = 0, theta2 = 0), DLL = "2d")

box_lower <- -4.5
box_upper <- 9.5
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
  labs(subtitle = "A: GHQ", size = "")

#' Shift by the mode
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = diag(c(1, 1)), dec.type = 2)

figA2 <- add_points(figA0, gg2) +
  labs(subtitle = "B: Shifted", size = "")

#' Adapt by the spectral
gg3 <- gg
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 1)

figA3 <- add_points(figA0, gg3) +
  labs(subtitle = "C: AGHQ", size = "")

#' PCA-AGHQ
gg4 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg4, m = mu, C = cov, dec.type = 1)

lambda <- eigen(cov)$values
cumsum(lambda) / sum(lambda)

xstart <- 6.1
ystart <- -2.3

x1end <- xstart + 4 * eigen(cov)$vectors[1, 1]
y1end <- ystart + 4 * eigen(cov)$vectors[2, 1]

x2end <- xstart + 1 * eigen(cov)$vectors[1, 2]
y2end <- ystart + 1 * eigen(cov)$vectors[2, 2]

figA4 <- add_points(figA0, gg4) +
  geom_segment(aes(x = xstart, y = ystart, xend = x1end, yend = y1end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x1end + 1, y = y1end - 3, label = "95%", col = "darkgrey") +
  geom_segment(aes(x = xstart, y = ystart, xend = x2end, yend = y2end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x2end, y = y2end - 2, label = "5%", col = "darkgrey") +
  labs(subtitle = "D: PCA-AGHQ", size = "")

ggsave(
  "fig2.png",
  plot = cowplot::plot_grid(figA1, figA2, figA3, figA4, ncol = 2),
  width = 6,
  height = 6,
  units = "in"
)

#' Figure 3

tmb <- readRDS("depends/tmb.rds")
outputs <- tmb$outputs

#' Adapted from figure from Naomi paper (Eaton 2021)
#' https://github.com/jeffeaton/naomi-model-paper/blob/master/analysis.R

indicators <- naomi:::add_output_labels(outputs) %>%
  left_join(
    select(outputs$meta_area, area_level, area_id, center_x, center_y),
    by = c("area_level", "area_id")
  ) %>%
  sf::st_as_sf() %>%
  mutate(
    area_level_label = fct_reorder(area_level_label, area_level),
    age_group_label = fct_reorder(age_group_label, as.integer(factor(age_group)))
  )

fig3adata <- indicators %>%
  filter(
    indicator == "prevalence",
    age_group == "Y015_049",
    sex == "both",
    calendar_quarter == "CY2016Q1",
    area_level == 4
  )

fig3a <- fig3adata %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  scale_fill_viridis_c(
    option = "C", direction = -1,
    begin = 0.1, end = 0.9,
    labels = scales::label_percent(1),
    breaks = c(0, 0.08, 0.16)
  ) +
  expand_limits(fill = 0) +
  labs(title = "A: HIV prevalence", fill = "") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8)),
    legend.text = element_text(size = rel(0.8)),
    legend.position = "bottom",
    legend.key.width = unit(1, "line")
  )

fig3bdata <- indicators %>%
  filter(
    indicator == "art_coverage",
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    sex == "both",
    area_level == 4
  )

fig3b <- fig3bdata %>%
  ggplot(aes(fill = mean)) +
  geom_sf(size = 0.1) +
  scale_fill_viridis_c(
    option = "D", direction = -1, begin = 0.1, end = 0.9,
    labels = scales::label_percent(1), limits = c(0.6, 0.853),
    breaks = c(0.6, 0.7, 0.8)
  ) +
  labs(title = "B: ART coverage", fill = "") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8)),
    legend.text = element_text(size = rel(0.8)),
    legend.position = "bottom",
    legend.key.width = unit(1, "line")
  )

fig3cdata <- indicators %>%
  filter(
    indicator %in% c("infections", "incidence"),
    age_group == "Y015_049",
    calendar_quarter == "CY2016Q1",
    sex == "both",
    area_level == 4
  )

fig3c <- fig3cdata %>%
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
  scale_size_area(max_size = 8, labels = scales::unit_format(unit = "K", scale = 1e-3, sep = "")) +
  labs(
    title = "C: HIV incidence",
    col = "",
    size = "",
    x = element_blank(),
    y = element_blank()
  ) +
  expand_limits(color = 0) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.8)),
    legend.text = element_text(size = rel(0.8)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(1, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-10, -10, -10, -10)
  )

ggsave(
  "fig3.png",
  plot = fig3a + fig3b + fig3c,
  width = 6,
  height = 5,
  units = "in"
)

#' Figure 4

#' Figure 5

exceedance <- readr::read_csv("depends/exceedance.csv")
exceedance_summary <- readr::read_csv("depends/exceedance-summary.csv")

fct_reorg <- function(fac, ...) {
  fct_recode(fct_relevel(fac, ...), ...)
}

fig5 <- exceedance %>%
  filter(indicator == "Second 90") %>%
  pivot_longer(cols = c("TMB", "aghq"), names_to = "method", values_to = "estimate") %>%
  mutate(
    method = fct_reorg(method, "TMB" = "TMB", "PCA-AGHQ" = "aghq")
  ) %>%
  ggplot(aes(x = tmbstan, y = estimate)) +
  geom_point(aes(col = sex), shape = 1, alpha = 0.6) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(values = c("#56B4E9", "#999999"), labels = c("Female", "Male")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(
    data = exceedance_summary %>%
      filter(indicator == "Second 90") %>%
      mutate(method = fct_relevel(method, "TMB", "PCA-AGHQ")),
    aes(x = -Inf, y = Inf, label = label),
    size = 3, hjust = 0, vjust = 1.5
  ) +
  facet_grid(indicator ~ method) +
  labs(x = "NUTS", y = "", col = "Sex", subtitle = "Estimated probabilities of meeting second 90 target") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  "fig5.png",
  plot = fig5,
  width = 6,
  height = 4,
  units = "in"
)
