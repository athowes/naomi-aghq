#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_sim")
# setwd("src/prev-anc-art_sim")

N_sim <- 1000

#' simplify = FALSE because I want a list of lists, not an array
sim_data <- replicate(N_sim, simulate_prev_anc_art(), simplify = FALSE)

saveRDS(sim_data, "sim_data.rds")

#' Add versions of data with different values of m
n <- 36
m_settings <- c(1, 10, 100)
for(m in m_settings){
  saveRDS(replicate(N_sim, simulate_prev_anc_art(m_prev = rep(m, n)), simplify = FALSE), paste0("sim_data_m", m, ".rds"))
}
