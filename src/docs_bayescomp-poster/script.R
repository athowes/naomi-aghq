# orderly::orderly_develop_start("docs_bayescomp-poster")
# setwd("src/docs_bayescomp-poster")

#' To-do
#' * [ ] Add MRC GIDA logo
#' * [ ] Add Oxford logo
#' * [ ] Add algorithm flowchart
#' * [ ]

rmarkdown::render("bayescomp-poster.Rmd")
pagedown::chrome_print("bayescomp-poster.html")
