#' Exclude parameters and data not required by naomi_simple model
local_exclude_inputs <- function(tmb_inputs) {

  simple_include_data <- c(
    "population_t1", "X_rho", "X_alpha", "X_lambda", "X_ancrho", "X_ancalpha",
    "Z_x", "Z_rho_x", "Z_rho_xs", "Z_rho_a", "Z_rho_as", "Z_rho_xa", "Z_alpha_x",
    "Z_alpha_xs", "Z_alpha_a", "Z_alpha_as", "Z_alpha_xa", "Z_lambda_x",
    "Z_ancrho_x", "Z_ancalpha_x", "log_asfr_t1_offset", "logit_anc_rho_t1_offset",
    "logit_anc_alpha_t1_offset", "logit_rho_offset", "logit_alpha_offset", "Q_x",
    "Q_x_rankdef", "n_nb", "adj_i", "adj_j", "Xgamma", "log_gamma_offset", "Xart_idx",
    "Xart_gamma", "omega", "OmegaT0", "sigma_OmegaT", "betaT0", "sigma_betaT",
    "ritaT", "X_15to49", "log_lambda_t1_offset", "X_15to49f", "X_paed_rho_ratio",
    "paed_rho_ratio_offset", "X_paed_lambda_ratio_t1", "x_prev", "n_prev",
    "A_prev", "x_artcov", "n_artcov", "A_artcov", "x_recent", "n_recent", "A_recent",
    "x_anc_prev_t1", "n_anc_prev_t1", "A_anc_prev_t1", "x_anc_artcov_t1", "n_anc_artcov_t1",
    "A_anc_artcov_t1", "A_artattend_t1", "x_artnum_t1", "A_artattend_mf",
    "A_art_reside_attend", "A_out", "A_anc_out", "calc_outputs", "report_likelihood"
  )

  simple_include_par <- c(
    "beta_rho", "beta_alpha", "beta_lambda", "beta_anc_rho", "beta_anc_alpha",
    "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs", "u_rho_a", "u_rho_as",
    "u_rho_xa", "ui_anc_rho_x", "ui_anc_alpha_x", "u_alpha_x", "us_alpha_x",
    "u_alpha_xs", "us_alpha_xs", "u_alpha_a", "u_alpha_as", "u_alpha_xa",
    "log_sigma_lambda_x", "ui_lambda_x", "logit_phi_rho_a", "log_sigma_rho_a",
    "logit_phi_rho_as", "log_sigma_rho_as", "logit_phi_rho_x", "log_sigma_rho_x",
    "logit_phi_rho_xs", "log_sigma_rho_xs", "log_sigma_rho_xa", "logit_phi_alpha_a",
    "log_sigma_alpha_a", "logit_phi_alpha_as", "log_sigma_alpha_as",
    "logit_phi_alpha_x", "log_sigma_alpha_x", "logit_phi_alpha_xs",
    "log_sigma_alpha_xs", "log_sigma_alpha_xa", "OmegaT_raw", "log_betaT",
    "log_sigma_ancrho_x", "log_sigma_ancalpha_x", "log_or_gamma", "log_sigma_or_gamma"
  )

  tmb_inputs$data <- tmb_inputs$data[simple_include_data]
  tmb_inputs$par_init <- tmb_inputs$par_init[simple_include_par]

  class(tmb_inputs) <- "naomi_simple_tmb_input"

  tmb_inputs
}

#' Create a quadrature grid
#'
#' @param dim The dimension of the grid
#' @param level The possible numbers of grid points per dimension
#' @param cut_off The cut-offs for standard deviations
#' @param sd A vector of length `dim` of standard deviations
sd_levels_ghe_grid <- function(dim, level, cut_off, sd) {
  stopifnot(length(level) == length(cut_off))
  stopifnot(dim == length(sd))
  levels <- vector(mode = "numeric", length = dim)
  for(i in seq_along(cut_off)) levels[sd > cut_off[i]] <- level[i]
  grid <- mvQuad::createNIGrid(dim = dim, "GHe", level = levels)
  grid
}
