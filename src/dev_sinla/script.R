#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

#' Create notebook
rmarkdown::render("sinla.Rmd")

#' Run stand-alone script
source("sinla.R")
