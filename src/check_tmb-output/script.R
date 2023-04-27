#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("check_tmb-output")
# setwd("src/check_tmb-output")

tmb <- readRDS("depends/tmb.rds")

qq_hist <- function(par, i) {
  df <- data.frame(y = tmb$fit$sample[[par]][i, ])

  qq <- ggplot(df, aes(sample = y)) +
    stat_qq(alpha = 0.2) +
    stat_qq_line() +
    labs(
      title = paste0("QQ plot for TMB ", par, "[", i, "]"),
      x = "", y = ""
    ) +
    theme_minimal()

  hist <- ggplot(df, aes(x = y)) +
    geom_histogram(alpha = 0.8) +
    labs(x = "", y = "") +
    theme_minimal()

  qq + hist
}

pdf("qq-checks.pdf", h = 4, w = 6.25)

qq_hist("ui_anc_alpha_x", 18)
qq_hist("ui_anc_alpha_x", 29)

dev.off()
