#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("multi-risk_sim")
# setwd("src/multi-risk_sim")

N_sim <- 1000

#' simplify = FALSE because I want a list of lists, not an array
sim_data <- replicate(N_sim, simulate_multi_risk(), simplify = FALSE)

saveRDS(sim_data, "sim_data.rds")
