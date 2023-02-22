library(patchwork)

#' Test code for inverting Hessian here

quad <- readRDS("depends/aghq.rds")$quad

#' Fit with k = 1 so there is only one modesandhessians entry
nrow(quad$modesandhessians)

H <- quad$modesandhessians$H[[1]]

H_df <- reshape2::melt(as.matrix(H))

pdf("H-matrix.pdf", h = 7, w = 6.25)

histogram <- ggplot(H_df, aes(x = log10(value))) +
  geom_histogram(alpha = 0.8) +
  labs(x = "log10(H[i, j])", y = "Count") +
  theme_minimal()

heatmap <- H_df %>%
  mutate(value_na = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = log10(value))) +
  labs(x = "", y = "", fill = "log10(H[i, j])") +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "i", y = "j") +
  theme_minimal()

heatmap / histogram

dev.off()
