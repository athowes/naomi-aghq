#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_process-results")
# setwd("src/prev-anc-art_process-results")

model0 <- readRDS("depends/results_model0.rds")
model1 <- readRDS("depends/results_model1.rds")
model2 <- readRDS("depends/results_model2.rds")
model3 <- readRDS("depends/results_model3.rds")
model4 <- readRDS("depends/results_model4.rds")

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

pdf("boxplots.pdf", h = 11, w = 8.5)

draw_boxplots(model0, c("tmb", "tmbstan")) + labs(title = "Model 0")
draw_boxplots(model1, c("tmb", "tmbstan", "aghq")) + labs(title = "Model 1")
draw_boxplots(model2, c("tmb", "tmbstan", "aghq")) + labs(title = "Model 2")
draw_boxplots(model3, c("tmb", "tmbstan", "aghq")) + labs(title = "Model 3")
draw_boxplots(model4, c("tmb", "tmbstan", "aghq")) + labs(title = "Model 4")

dev.off()
