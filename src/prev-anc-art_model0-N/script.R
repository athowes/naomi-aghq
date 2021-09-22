#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model0-N")
# setwd("src/prev-anc-art_model0-N")

sim_data <- readRDS("depends/sim_data.rds")

#' TMB preparation
compile("model0.cpp")
dyn.load(dynlib("model0"))

#' Run for each simulated dataset (using head for now for speed)
results <- map(head(sim_data), run_model0)

#' Save to artefact
saveRDS(results, "results.rds")
