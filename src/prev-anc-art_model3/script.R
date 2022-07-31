#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model3")
# setwd("src/prev-anc-art_model3")

sim_data <- readRDS("depends/sim_data.rds")

#' TMB preparation
compile("model3.cpp")
dyn.load(dynlib("model3"))

#' Run for each simulated dataset (using head for now for speed)
results <- map(head(sim_data), run_model3)

#' Save to artefact
saveRDS(results, "results.rds")
