# orderly::orderly_develop_start("docs_paper-jtb")
# setwd("src/docs_paper-jtb")

source("figs.R")

#' Conversion of figures from .pdf at specified resolution
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/homebrew/bin", sep = ":"))

convert_pdf_png <- function(name) {
  command <- paste0(
    "convert -density 300 depends/", name, ".pdf -scene 1 -background white",
    " -alpha remove -alpha off -quality 80 depends/", name, ".png"
  )
  system(command)
}

#' For paper.Rmd
convert_pdf_png("nodes-samples-comparison")

#' Render documents
rmarkdown::render("paper.Rmd", clean = FALSE)
rmarkdown::render("appendix.Rmd", clean = FALSE)

pdftools::pdf_combine(c("paper.pdf", "appendix.pdf"), output = "naomi-aghq.pdf")

#' I just used this function to check some things as I was writing!
kish_ess <- function(w) {
  sum(w)^2 / sum(w^2)
}
