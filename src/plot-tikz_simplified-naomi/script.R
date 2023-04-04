#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("plot-tikz_simplified-naomi")
# setwd("src/plot-tikz_simplified-naomi")

tools::texi2dvi("simplified-naomi.tex", pdf = TRUE, clean = TRUE)

convert_pdf_png <- function(name) {
  command <- paste0(
    "convert -density 300 ", name, ".pdf -scene 1 -background white",
    " -alpha remove -alpha off -quality 80 ", name, ".png"
  )
  system(command)
}

convert_pdf_png("simplified-naomi")
