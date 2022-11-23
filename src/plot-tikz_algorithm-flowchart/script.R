#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("plot-tikz_algorithm-flowchart")
# setwd("src/plot-tikz_algorithm-flowchart")

tools::texi2dvi("algorithm-flowchart.tex", pdf = TRUE, clean = TRUE)
