#' Figure 1

#' A screenshot from https://naomi.unaids.org/

#' Figure 2

mu <- c(1, 1.5)
cov <- matrix(c(2, 1, 1, 1), ncol = 2)

obj <- function(theta) {
  mvtnorm::dmvnorm(theta, mean = mu, sigma = cov)
}

grid <- expand.grid(
  theta1 = seq(-2, 5, length.out = 700),
  theta2 = seq(-2, 5, length.out = 700)
)

ground_truth <- cbind(grid, pdf = obj(grid))

fig2base <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = "lightgrey") +
  coord_fixed(xlim = c(-2, 4.5), ylim = c(-2, 4.5), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal() +
  guides(size = FALSE) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = grid::unit(c(0,0,0,0), "mm")
  )

gg <- mvQuad::createNIGrid(2, "GHe", 3)

add_points <- function(fig2base, gg) {
  fig2base +
    geom_point(
      data = mvQuad::getNodes(gg) %>%
        as.data.frame() %>%
        mutate(weights = mvQuad::getWeights(gg)),
      aes(x = V1, y = V2, size = weights),
      alpha = 0.8,
      col = "#009E73",
      inherit.aes = FALSE
    ) +
    scale_size_continuous(range = c(1, 2))
}

fig2a <- add_points(fig2base, gg) +
  labs(title = "GHQ")

#' Adapt by the mean
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = diag(c(1, 1)), dec.type = 2)

fig2b <- add_points(fig2base, gg2) +
  labs(title = "Shift")

#' Adapt by the spectral
gg3 <- gg
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 1)

fig2c <- add_points(fig2base, gg3) +
  labs(title = "Rotate")

#' PCA-AGHQ
gg4 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg4, m = mu, C = cov, dec.type = 1)

fig2d <- add_points(fig2base, gg4) +
  labs(title = "Keep")

ggsave(
  "fig2.png",
  plot = (fig2a + fig2b) / (fig2c + fig2d),
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
    labels = scales::label_percent(1)
  ) +
  expand_limits(fill = 0) +
  labs(title = "HIV prevalence", fill = "") +
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
    legend.key.width = unit(1.5, "line")
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
    option = "D", direction = -1, begin = 0.05, end = 0.9,
    labels = scales::label_percent(1), limits = c(0.6, 0.853)
  ) +
  labs(title = "ART coverage", fill = "") +
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
    legend.key.width = unit(1.5, "line")
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
  scale_size_area(max_size = 8) +
  labs(
    title = "HIV incidence",
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
    legend.key.width = unit(1.5, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-10, -10, -10, -10)
  )

ggsave(
  "fig3.png",
  plot = fig3a + fig3b + fig3c,
  width = 6,
  height = 6,
  units = "in"
)
