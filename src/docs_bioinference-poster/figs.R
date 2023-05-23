#' Figure 1

#' A screenshot from https://naomi.unaids.org/

#' Figure 2

#' Created in the poster

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
  labs(fill = "HIV prevalence") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
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
  labs(fill = "ART coverage") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
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
    color = "Incidence per 1000\nperson years",
    size = "Annual new infections",
    x = element_blank(),
    y = element_blank()
  ) +
  expand_limits(color = 0) +
  coord_sf(expand = FALSE) +
  theme_minimal(10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    strip.text = element_text(face = "bold", size = rel(0.85)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)),
    plot.tag = element_text(face = "bold")
  )

ggsave(
  "fig3.png",
  plot = fig3a + fig3b + fig3c,
  width = 11,
  height = 6,
  units = "in"
)
