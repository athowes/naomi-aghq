# orderly::orderly_develop_start("docs_bayescomp-poster")
# setwd("src/docs_bayescomp-poster")

#' Create figures
source("figs.R")

#' Create poster
rmarkdown::render("bayescomp-poster.Rmd")
pagedown::chrome_print("bayescomp-poster.html")
