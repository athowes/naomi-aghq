#' A local version of naomi::make_tmb_obj, edited to work with DLL = "naomi_simple"
local_make_tmb_obj <- function(data, par, calc_outputs = 0L, inner_verbose, progress = NULL, map = NULL, DLL = "naomi_simple") {
  # Begin expose naomi:::make_tmb_obj
  # https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L496
  data$calc_outputs <- as.integer(calc_outputs)

  integrate_out <- c(
    "beta_rho", "beta_alpha", "beta_alpha_t2", "beta_lambda", "beta_asfr",
    "beta_anc_rho", "beta_anc_alpha", "beta_anc_rho_t2", "beta_anc_alpha_t2",
    "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs", "u_rho_a", "u_rho_as",
    "u_rho_xa", "u_alpha_x", "us_alpha_x", "u_alpha_xs", "us_alpha_xs",
    "u_alpha_a", "u_alpha_as", "u_alpha_xt", "u_alpha_xa", "u_alpha_xat",
    "u_alpha_xst", "ui_lambda_x", "logit_nu_raw", "ui_asfr_x", "ui_anc_rho_x",
    "ui_anc_alpha_x", "ui_anc_rho_xt", "ui_anc_alpha_xt", "log_or_gamma",
    "log_or_gamma_t1t2"
  )

  if (DLL == "naomi_simple") {
    exclude_random_pars <- c(
      "beta_alpha_t2", "beta_asfr", "beta_anc_rho_t2", "beta_anc_alpha_t2",
      "u_alpha_xt", "u_alpha_xat", "u_alpha_xst", "logit_nu_raw", "ui_asfr_x",
      "ui_anc_rho_xt", "ui_anc_alpha_xt", "log_or_gamma_t1t2"
    )

    integrate_out <- setdiff(integrate_out, exclude_random_pars)
  }

  if(DLL == "naomi_simple_x_index") {
    integrate_out <- "x_minus_i"
  }

  obj <- TMB::MakeADFun(
    data = data,
    parameters = par,
    DLL = DLL,
    silent = !inner_verbose,
    random = integrate_out,
    map = map
  )

  if (!is.null(progress)) {
    obj$fn <- naomi:::report_progress(obj$fn, progress)
  }

  obj
  # End expose naomi:::make_tmb_obj
}

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

#' Inference for the Naomi model using aghq, edited to work with DLL = "naomi_simple"
fit_aghq <- function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, DLL = "naomi_simple", ...) {
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL = DLL)
  quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, ...)
  objout <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 1L, inner_verbose, progress, map, DLL = DLL)
  quad$obj <- objout
  quad
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

#' Create a function to do the PCA rescaling, which also adapts according to the mean and reweights the nodes:
#'
#' @param m Mean vector
#' @param C Covariance matrix
#' @param s Small grid dimension
#' @param k Number of points per small grid dimension
pca_rescale <- function(m, C, s, k) {
  d <- nrow(C)
  stopifnot(d == length(m))
  eigenC <- eigen(C)
  lambda <- eigenC$values
  E <- eigenC$vectors
  E_s <- E[, 1:s]
  gg_s <- mvQuad::createNIGrid(dim = s, type = "GHe", level = k)
  nodes_out <- t(E_s %*% diag(sqrt(lambda[1:s]), ncol = s) %*% t(mvQuad::getNodes(gg_s)))
  for(j in 1:d) nodes_out[, j] <- nodes_out[, j] + m[j] # Adaption
  weights_out <- mvQuad::getWeights(gg_s) * as.numeric(mvQuad::getWeights(mvQuad::createNIGrid(dim = d - s, type = "GHe", level = 1)))
  weights_out <- det(chol(C)) * weights_out # Adaption

  # Putting things into a mvQuad format manually
  gg <- mvQuad::createNIGrid(dim = d, type = "GHe", level = 1)
  gg$level <- rep(NA, times = d)
  gg$ndConstruction <- "PCA"
  gg$nodes <- nodes_out
  gg$weights <- weights_out
  return(gg)
}

plot_matrix <- function(M) {
  name <- deparse(substitute(M))
  M_df <- reshape2::melt(as.matrix(M))

  ggplot(M_df, aes(x = Var1, y = Var2, fill = value)) +
    labs(x = "", y = "", fill = paste0(name, "[i, j]")) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x = "i", y = "j") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

plot_total_variation <- function(eigen, label_x) {
  ggplot(data = NULL, aes(x = 1:length(eigen$values), y = cumsum(eigen$values) / sum(eigen$values))) +
    geom_point() +
    geom_hline(yintercept = 0.9, col = "grey", linetype = "dashed") +
    annotate("text", x = label_x, y = 0.875, label = "90% of total correlation explained", col = "grey") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "PCA dimensions included", y = "Total variation explained") +
    theme_minimal()
}

plot_pc_loadings <- function(eigen) {
  reshape2::melt(as.matrix(eigen$vectors)) %>%
    ggplot(aes(x = Var1, y = factor(Var2), fill = value)) +
    labs(x = "", y = "", fill = "") +
    geom_tile() +
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#E69F00", "white", "#009E73"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0),
      limits = c(-1.1, 1.1)
    ) +
    labs(x = "Hyper", y = "Principal component loading") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
