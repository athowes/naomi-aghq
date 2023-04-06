#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_scale-grid")
# setwd("src/dev_scale-grid")

#' Create notebook
rmarkdown::render("scale-grid.Rmd")

#' Run astronomy example
rmarkdown::render("astro.Rmd")
