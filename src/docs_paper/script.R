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
#' Found a way to cross-link references: https://stackoverflow.com/questions/52531637/knitr-rmarkdown-latex-how-to-cross-reference-figures-and-tables-in-2-different/52532269#52532269
rmarkdown::render("paper.Rmd", clean = FALSE)
rmarkdown::render("appendix.Rmd", clean = FALSE)

tinytex::pdflatex("paper.tex", clean = FALSE)
tinytex::pdflatex("appendix.tex", clean = FALSE)

tinytex::pdflatex("paper.tex")
tinytex::pdflatex("appendix.tex")

#' I just used this function to check some things as I was writing!
kish_ess <- function(w) {
  sum(w)^2 / sum(w^2)
}
