simulate_prev_anc_art <- function(n = 30,
                                  m_prev = rep(250, n),
                                  beta_prev = -2.4,
                                  tau_phi_prev = 4, #' sigma_phi_prev is 0.5
                                  m_anc = rep(10^4, n),
                                  beta_anc = -0.2,
                                  tau_b_anc = 100, #' sigma_b_anc is 0.1
                                  N_art = rep(10^5, n),
                                  beta_art = 0.7,
                                  tau_phi_art = (1 / 0.35)^2) { #' sigma_phi_art is 0.35

  #' Survey prevalence
  eta_prev <- beta_prev + rnorm(n, 0, 1 / sqrt(tau_phi_prev))
  rho_prev <- plogis(eta_prev)
  y_prev <- rbinom(n, m_prev, rho_prev)

  #' ANC
  b_anc <- beta_anc + rnorm(n, 0,  1 / sqrt(tau_b_anc))
  eta_anc <- eta_prev + b_anc
  rho_anc <- plogis(eta_anc)
  y_anc <- rbinom(n, m_anc, rho_anc)

  #' ART
  eta_art <- beta_art + rnorm(n, 0, 1 / sqrt(tau_phi_art))
  alpha_art <- plogis(eta_art)
  A_art <- rbinom(n, floor(N_art * rho_prev), alpha_art)

  return(list(
    n = n,
    m_prev = m_prev,
    m_anc = m_anc,
    N_art = N_art,
    eta_prev = eta_prev,
    rho_prev = rho_prev,
    y_prev = y_prev,
    b_anc = b_anc,
    eta_anc = eta_anc,
    rho_anc = rho_anc,
    y_anc = y_anc,
    eta_art = eta_art,
    alpha_art = alpha_art,
    A_art = A_art
  ))
}
