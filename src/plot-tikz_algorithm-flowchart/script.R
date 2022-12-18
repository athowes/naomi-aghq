#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("plot-tikz_algorithm-flowchart")
# setwd("src/plot-tikz_algorithm-flowchart")

tools::texi2dvi("algorithm-flowchart.tex", pdf = TRUE, clean = TRUE)

system("convert -density 300 -background white -alpha background -alpha off algorithm-flowchart.pdf algorithm-flowchart.png")

