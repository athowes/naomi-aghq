#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("docs_11-04-23-explainer")
# setwd("src/docs_11-04-23-explainer")

cbpalette <- multi.utils::cbpalette()

rmarkdown::render("11-04-23-explainer.Rmd")
