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

#' A local version of naomi::fit_tmb, edited to work with DLL = "naomi_simple"
local_fit_tmb <- function(tmb_input, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL, map = NULL, DLL = "naomi") {
  # Begin expose naomi::fit_tmb
  # https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L557

  if (DLL == "naomi_simple") {
    stopifnot(inherits(tmb_input, "naomi_simple_tmb_input"))
  } else {
    stopifnot(inherits(tmb_input, "naomi_tmb_input"))
  }

  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL = DLL)
  data <- tmb_input$data
  par <- tmb_input$par_init
  calc_outputs <- 0L

  trace <- ifelse(outer_verbose, 1, 0)

  f <- withCallingHandlers(
    stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace, iter.max = max_iter)),
    warning = function(w) {
      if (grepl("NA/NaN function evaluation", w$message)) invokeRestart("muffleWarning")
    }
  )

  if (f$convergence != 0) warning(paste("convergence error:", f$message))

  if (outer_verbose) message(paste("converged:", f$message))

  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  objout <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 1L, inner_verbose, DLL = DLL)
  f$mode <- objout$report(f$par.full)
  val <- c(f, obj = list(objout))
  class(val) <- "naomi_fit"

  # Returns parameter estimate (mode) and Hessian
  val

  # End expose naomi::fit_tmb
}

#' A local version of naomi::sample_tmb
local_sample_tmb <- function(fit, nsample = 1000, rng_seed = NULL, random_only = TRUE, verbose = FALSE) {
  # Begin expose naomi::sample_tmb
  # https://github.com/mrc-ide/naomi/blob/65ac94517b910ac517a45f41e824824e1907a3c4/R/tmb-model.R#L624
  set.seed(rng_seed)
  stopifnot(methods::is(fit, "naomi_fit"))
  stopifnot(nsample > 1)
  to_tape <- TMB:::isNullPointer(fit$obj$env$ADFun$ptr)
  if (to_tape) fit$obj$retape(FALSE)
  if (!random_only) {
    if (verbose) print("Calculating joint precision")
    hess <- naomi:::sdreport_joint_precision(fit$obj, fit$par.fixed)
    if (verbose) print("Inverting precision for joint covariance")
    cov <- solve(hess)
    if (!isSymmetric(cov, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
      stop("cov must be a symmetric matrix")
    }
    if (verbose) print("Drawing sample")
    smp <- mvtnorm::rmvnorm(n = nsample, mean = fit$par.full, sigma = cov, method = "eigen", checkSymmetry = FALSE)
  } else {
    r <- fit$obj$env$random # Indices of the random effects
    par_f <- fit$par.full[-r] # Mode of the fixed effects
    par_r <- fit$par.full[r] # Mode of the random effects
    hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE) # Hessian of the random effects
    smp_r <- naomi:::rmvnorm_sparseprec(nsample, par_r, hess_r) # Sample from the random effects
    smp <- matrix(0, nsample, length(fit$par.full)) # Create data structure for sample storage
    smp[, r] <- smp_r # Store random effect samples
    smp[, -r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE) # For fixed effects store mode
    colnames(smp)[r] <- colnames(smp_r) # Random effect names
    colnames(smp)[-r] <- names(par_f) # Fixed effect names
  }
  if (verbose) print("Simulating outputs")

  # Number of rows nrow(smp) equal number of samples
  # Number of columns ncol(smp) equal number of parameters (latent field and hyper)
  # sum(lengths(tmb_input$par_init)) is the same as ncol(smp)

  # Given samples, use TMB template to generate output samples
  sim <- apply(smp, 1, fit$obj$report)
  r <- fit$obj$report()
  if (verbose) print("Returning sample")
  fit$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
  is_vector <- vapply(fit$sample, inherits, logical(1), "numeric")
  fit$sample[is_vector] <- lapply(fit$sample[is_vector], matrix, nrow = 1)
  names(fit$sample) <- names(r)

  fit
  # End expose naomi::sample_tmb
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
