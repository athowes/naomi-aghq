# orderly::orderly_develop_start("docs_statml-poster")
# setwd("src/docs_statml-poster")

#' Create figures
source("figs.R")

#' Create poster
rmarkdown::render("statml-poster.Rmd")
pagedown::chrome_print("statml-poster.html")
