#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

#' Create notebook
rmarkdown::render("sinla.Rmd")

cbpalette <- multi.utils::cbpalette()

sim_data_m1 <- readRDS("depends/sim_data_m1.rds")
sim_data_m10 <- readRDS("depends/sim_data_m10.rds")
sim_data_m100 <- readRDS("depends/sim_data_m100.rds")
sim_data_m250 <- readRDS("depends/sim_data.rds")

run_analysis(sim_data_m250, "m250")
