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

#' Figure 4
fig4data <- readRDS("depends/beta_anc_rho.rds") %>%
  filter(method != "aghq")

colours <- c("#56B4E9", "#009E73", "#E69F00")

mean <- fig4data %>%
  filter(method == "tmbstan") %>%
  summarise(mean = mean(samples)) %>%
  pull(mean) %>%
  round(digits = 3)

sd <- fig4data %>%
  filter(method == "tmbstan") %>%
  summarise(sd = sd(samples)) %>%
  pull(sd) %>%
  round(digits = 3)

fig4a <- ggplot(fig4data, aes(x = samples, fill = method, col = method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  facet_grid(method~.) +
  labs(x = "beta_anc_rho", y = "Density", fill = "Method") +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none")

grid <- seq(from = min(fig4data$samples), to = max(fig4data$samples), length.out = 1000)

tmb_ecdf <- stats::ecdf(filter(fig4data, method == "TMB") %>% pull(samples))
tmb_ecdf_df <- data.frame(x = grid, ecdf = tmb_ecdf(grid), method = "TMB")

adam_ecdf <- stats::ecdf(filter(fig4data, method == "adam") %>% pull(samples))
adam_ecdf_df <- data.frame(x = grid, ecdf = adam_ecdf(grid), method = "adam")

tmbstan_ecdf <- stats::ecdf(filter(fig4data, method == "tmbstan") %>% pull(samples))
tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "tmbstan")

# Add ECDF differences
tmb_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - tmb_ecdf_df$ecdf
adam_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - adam_ecdf_df$ecdf
tmbstan_ecdf_df$ecdf_diff <- 0

absmax <- function(x) x[which.max(abs(x))]

ks_tmb <- absmax(tmb_ecdf_df$ecdf_diff)
ks_adam <- absmax(adam_ecdf_df$ecdf_diff)

ecdf_df <- bind_rows(tmb_ecdf_df, adam_ecdf_df, tmbstan_ecdf_df)

ecdf_df$method <- factor(ecdf_df$method, levels = c("TMB", "adam", "tmbstan"))

ks_labeller <- function(x) toString(round(abs(x), 2))

fig4b <- ggplot(ecdf_df, aes(x = x, y = ecdf_diff, col = method)) +
  geom_line() +
  geom_abline(intercept = ks_tmb, slope = 0, col = colours[1], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_tmb, label = ks_labeller(ks_tmb), col = colours[1], alpha = 0.8) +
  geom_abline(intercept = ks_adam, slope = 0, col = colours[2], linetype = "dashed", alpha = 0.8) +
  annotate("text", x = 1.1 * max(ecdf_df$x), y = ks_adam, label = ks_labeller(ks_adam), col = colours[2], alpha = 0.8) +
  scale_color_manual(values = colours) +
  labs(x = "beta_anc_rho", y = "ECDF difference from tmbstan") +
  guides(col = "none") +
  coord_cartesian(xlim = c(min(ecdf_df$x), max(ecdf_df$x)), clip = "off") +
  theme_minimal() +
  theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))

ggsave(
  "fig4.png",
  plot = fig4a + fig4b,
  width = 11,
  height = 6,
  units = "in"
)
