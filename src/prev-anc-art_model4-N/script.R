#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model4-N")
# setwd("src/prev-anc-art_model4-N")

sim_data <- readRDS("depends/sim_data.rds")

#' TMB preparation
compile("model4.cpp")
dyn.load(dynlib("model4"))

#' Run for each simulated dataset (using head for now for speed)
results <- map(head(sim_data), run_model4)

#' Save to artefact
saveRDS(results, "results.rds")
