#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model0")
# setwd("src/prev-anc-art_model0")

sim_data <- readRDS("depends/sim_data.rds")

#' TMB preparation
compile("model0.cpp")
dyn.load(dynlib("model0"))

#' TODO
#' sd_out has no jointPrecision because there are no random effects in this model
#' Need to extend the functionality of inf.utils::sample_tmb to account for this case
#' Try rerunning this once that is done!

#' Run for each simulated dataset (using head for now for speed)
results <- map(head(sim_data), run_model0)

#' Save to artefact
saveRDS(results, "results.rds")
