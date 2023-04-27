# orderly::orderly_develop_start("docs_bioinference-poster")
# setwd("src/docs_bioinference-poster")

#' Create figures
source("figs.R")

#' Create poster
rmarkdown::render("bioinference-poster.Rmd")
pagedown::chrome_print("bioinference-poster.html")
