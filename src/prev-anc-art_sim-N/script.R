#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_sim-N")
# setwd("src/prev-anc-art_sim-N")

N_sim <- 1000

#' simplify = FALSE because I want a list of lists, not an array
sim_data <- replicate(N_sim, simulate_prev_anc_art(), simplify = FALSE)

saveRDS(sim_data, "sim_data.rds")
