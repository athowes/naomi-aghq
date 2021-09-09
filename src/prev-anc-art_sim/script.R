#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_sim")
# setwd("src/prev-anc-art_sim")

set.seed(1)

#' Simulate data

#' Survey prevalence
n <- 5
m_prev <- rep(20, n)
beta_prev <- -1
tau_phi_prev <- 1
eta_prev <- beta_prev + rnorm(n, 0, 1 / sqrt(tau_phi_prev))
rho_prev <- plogis(eta_prev)
y_prev <- rbinom(n, m_prev, rho_prev)

#' ANC
m_anc <- rep(50, n)
beta_anc <- 0.5
tau_b_anc <- 4
b_anc <- beta_anc + rnorm(n, 0,  1 / sqrt(tau_b_anc))
eta_anc <- eta_prev + b_anc
rho_anc <- plogis(eta_anc)
y_anc <- rbinom(n, m_anc, rho_anc)

#' ART
N_art <- rep(500, n)
beta_art <- 0
tau_phi_art <- 1
eta_art <- beta_art + rnorm(n, 0, 1 / sqrt(tau_phi_art))
alpha_art <- plogis(eta_art)
A_art <- rbinom(n, floor(N_art * rho_prev), alpha_art)

#' Save image of workspace
save.image(file = "sim_data.rdata")
