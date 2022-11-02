#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("docs_15-11-22-seminar")
# setwd("src/docs_15-11-22-seminar")

cbpalette <- multi.utils::cbpalette()

rmarkdown::render("15-11-22-seminar.Rmd")
