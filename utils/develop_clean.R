#' Set the working directory to the project root
setwd(rprojroot::find_rstudio_root_file())

#' List of all the reports
tasks <- list.files("src")

#' Remove all artefacts, depends
lapply(tasks, orderly::orderly_develop_clean)
