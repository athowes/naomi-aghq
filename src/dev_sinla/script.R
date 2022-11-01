#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

cbpalette <- multi.utils::cbpalette()

#' Create notebook
rmarkdown::render("sinla.Rmd")

#' Experiments
sim_data_m1 <- readRDS("depends/sim_data_m1.rds")
sim_data_m10 <- readRDS("depends/sim_data_m10.rds")
sim_data_m100 <- readRDS("depends/sim_data_m100.rds")
sim_data_m250 <- readRDS("depends/sim_data.rds")

#' IID random effects
run_analysis_model1(sim_data_m1, "m1")
run_analysis_model1(sim_data_m10, "m10")
run_analysis_model1(sim_data_m100, "m100")
run_analysis_model1(sim_data_m250, "m250")

#' ICAR random effects
# run_analysis_model1_icar(sim_data_m1, "m1")
# run_analysis_model1_icar(sim_data_m10, "m10")
# run_analysis_model1_icar(sim_data_m100, "m100")
# run_analysis_model1_icar(sim_data_m250, "m250")
