# orderly::orderly_develop_start("docs_paper")
# setwd("src/docs_paper")

source("figures.R")

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
rmarkdown::render("paper.Rmd")
rmarkdown::render("appendix.Rmd")


kish_ess <- function(w) {
  sum(w)^2 / sum(w^2)
}

kish_ess(c(1, 1, 1, 1, 1))
kish_ess(c(2, 5, 2, 2, 2))
