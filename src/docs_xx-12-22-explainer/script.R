#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("docs_xx-12-22-explainer")
# setwd("src/docs_xx-12-22-explainer")

cbpalette <- multi.utils::cbpalette()

rmarkdown::render("xx-12-22-explainer.Rmd")
