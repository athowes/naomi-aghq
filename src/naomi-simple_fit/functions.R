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
  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)
  quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, ...)
  quad$obj <- obj
  quad
}

#' Uncertainty for the Naomi model using aghq
sample_aghq <- function(quad, M, verbose = TRUE) {
  # Note that with k = 1, sample_marginal just returns the mode for the hypers
  # This needs to be debugged, maybe?
  if (verbose) print("Sampling from aghq")
  samp <- aghq::sample_marginal(quad, M)

  # This part replaces samples from TMB with samples from aghq
  if (verbose) print("Rearranging samples")
  r <- quad$obj$env$random
  smp <- matrix(0, M, length(quad$obj$env$par))
  smp[, r] <- unname(t(samp$samps))
  names(samp$thetasamples) <- names(samp$theta)
  smp[, -r] <- unname(as.matrix(bind_rows(samp$thetasamples)))
  smp <- as.data.frame(smp)
  colnames(smp)[r] <- rownames(samp$samps)
  colnames(smp)[-r] <- names(samp$thetasamples)

  # This part is the same as TMB
  if (verbose) print("Simulating from model")
  sim <- apply(smp, 1, quad$obj$report)
  r <- quad$obj$report()
  if (verbose) print("Returning sample")
  quad$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
  is_vector <- vapply(quad$sample, inherits, logical(1), "numeric")
  quad$sample[is_vector] <- lapply(quad$sample[is_vector], matrix, nrow = 1)
  names(quad$sample) <- names(r)

  quad
}

#' Inference for the Naomi model using tmbstan, edited to work with DLL = "naomi_simple"
fit_tmbstan <- function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, DLL = "naomi_simple", ...) {
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)
  fit <- tmbstan::tmbstan(obj, ...)
  tmbstanfit <- setClass("tmbstanfit", contains = "stanfit", slots = c(obj = "list"))
  fit <- as(fit, "tmbstanfit")
  fit@obj <- obj
  fit
}

#' Uncertainty for the Naomi model using tmbstan
sample_tmbstan <- function(mcmc, M, verbose = TRUE) {
  if (verbose) print("Samples already available for tmbstan!")
  if (verbose) print("Simulating from model")
}

#' Inference for the Naomi model using aghq plus Laplace marginals, edited to work with DLL = "naomi_simple"
fit_adam <- function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, DLL = "naomi_simple",  ...) {
  warning("For the aghq, k must be set to 1 as thing stand currently.")
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)
  quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = 1)
  quad$obj <- obj

  random <- obj$env$random
  x_names <- names(obj$env$par[random])
  x_lengths <- lengths(tmb_input$par_init[unique(x_names)])
  x_starts <- cumsum(x_lengths) - x_lengths

  tmb_input_i <- tmb_input
  tmb_input_i$par_init$x_minus_i <- rep(0, sum(x_lengths) - 1)
  tmb_input_i$par_init$x_i <- 0
  tmb_input_i$par_init[unique(x_names)] <- NULL
  tmb_input_i$data$x_lengths <- x_lengths
  tmb_input_i$data$x_starts <- x_starts

  theta_names <- make.unique(names(obj$par), sep = "")
  theta <- as.numeric(quad$modesandhessians[theta_names])

  .f <- function(i) {
    tmb_input_i$data$i <- i
    obj_i <- local_make_tmb_obj(tmb_input_i$data, tmb_input_i$par, calc_outputs = 0L, inner_verbose, progress, DLL = "naomi_simple_x_index")
    random_i <- obj_i$env$random
    mode_i <- quad$modesandhessians[["mode"]][[1]][-i]


    gg <- create_approx_grid(quad$modesandhessians, i = i, k = 5)

    out <- data.frame(
      index = i,
      par = x_names[i],
      x = mvQuad::getNodes(gg),
      w = mvQuad::getWeights(gg)
    )

    # Inclusion of the weights argument here is meaningless for now, but will become
    # of importance when there are multiple hyperparameter grid points used. We also
    # always put mode_i in as the last parameter so that optimisation starts there.
    # It's unclear to me if there is a better way to do this.
    .g <- function(x) {
      obj_i$env$last.par[random_i] <- mode_i
      return(as.numeric(- obj_i$fn(c(x, theta))) + log(quad$normalized_posterior$nodesandweights$weights))
    }

    out$lp <- purrr::map_dbl(out$x, .g) # Calculate log-posterior
    lognormconst <- logSumExpWeights(out$lp, out$w) # Calculate log normalising constant
    out$lp_normalised <- out$lp - lognormconst # Calculate normalised log-posterior

    return(out)
  }

  d <- length(random)
  laplace_marginals <- purrr::map(.x = 1:d, .f = .f)
  laplace_marginals <- dplyr::bind_rows(laplace_marginals)

  out <- list("quad" = quad, "laplace_marginals" = laplace_marginals)
  return(out)
}

#' Uncertainty for the Naomi model using adam
sample_adam <- function(adam, M, verbose = TRUE) {
  if (verbose) print("Sampling from adam")
  r <- adam$quad$obj$env$random
  d <- length(r)
  x_names <- names(adam$quad$obj$env$par[r])
  theta_names <- names(adam$quad$obj$env$par[-r])
  samps <- matrix(data = NA, nrow = d, ncol = M)
  rownames(samps) <- x_names

  # Mirroring the aghq data structures as in aghq::sample_marginal.marginallaplace
  # Might want to change these as well as in the aghq package eventually

  # Laplace latent field marginals
  for(i in 1:d) {
    marginal <- adam$laplace_marginals[adam$laplace_marginals$index == i, ]
    # Extra normalisation check in here to avoid any issues with slightly not integrating to one
    pdf_and_cdf <- compute_pdf_and_cdf(marginal$x, marginal$lp_normalised, normalise = TRUE)
    samps[i, ] <- sample_cdf(pdf_and_cdf, M = M)
  }

  # Hyperparameter marginals
  thetasamples <- list()
  for(j in 1:length(adam$quad$marginals)) {
    marginal <- adam$quad$marginals[[j]]
    colnames(marginal)[grep("theta", colnames(marginal))] <- "theta"
    # Don't normalise because maybe just one point (for now) and the trapezoid
    # rule will break
    pdf_and_cdf <- compute_pdf_and_cdf(marginal$theta, marginal$logmargpost)
    thetasamples[[j]] <- unname(sample_cdf(pdf_and_cdf, M = M))
  }

  samp <- list()
  samp$samps <- samps
  samp$thetasamples <- thetasamples

  # This part replaces samples from TMB with samples from aghq
  if (verbose) print("Rearranging samples")
  smp <- matrix(0, M, length(adam$quad$obj$env$par))
  smp[, r] <- unname(t(samp$samps))
  names(samp$thetasamples) <- theta_names
  smp[, -r] <- unname(as.matrix(bind_rows(samp$thetasamples)))
  smp <- as.data.frame(smp)
  colnames(smp)[r] <- x_names
  colnames(smp)[-r] <- theta_names

  # This part is the same as TMB
  if (verbose) print("Simulating from model")
  sim <- apply(smp, 1, adam$quad$obj$report)
  r <- adam$quad$obj$report()
  if (verbose) print("Returning sample")
  adam$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r), numeric), names(r))
  is_vector <- vapply(adam$sample, inherits, logical(1), "numeric")
  adam$sample[is_vector] <- lapply(adam$sample[is_vector], matrix, nrow = 1)
  names(adam$sample) <- names(r)

  adam
}

#' Version of output_package() to extract only T1
local_output_package_naomi_simple <- function(naomi_fit, naomi_data, na.rm = FALSE) {

  stopifnot(is(naomi_fit, "naomi_fit"))
  stopifnot(is(naomi_data, "naomi_data"))

  indicators <- local_extract_indicators_naomi_simple(naomi_fit, naomi_data, na.rm = na.rm)

  art_attendance <- local_extract_art_attendance_naomi_simple(naomi_fit, naomi_data, na.rm = na.rm)

  meta_area <- naomi_data$areas %>%
    dplyr::filter(area_id %in% unique(naomi_data$mf_out$area_id)) %>%
    dplyr::select(area_level, area_level_label, area_id, area_name,
                  parent_area_id, spectrum_region_code, area_sort_order,
                  center_x, center_y, geometry) %>%
    sf::st_as_sf()

  meta_period <- naomi:::get_period_metadata(naomi_data$calendar_quarter1)
  meta_age_group <- naomi:::get_age_groups()

  # Fitting outputs
  fit <- list()
  fit$model_options <- naomi_data$model_options
  fit$data_options <- naomi_data$data_options
  fit$calibration_options <- naomi_data$calibration_options
  fit$spectrum_calibration <- naomi_data$spectrum_calibration

  meta_indicator <- get_meta_indicator()
  meta_indicator <- dplyr::filter(meta_indicator, indicator %in% indicators$indicator)

  inputs_outputs <- align_inputs_outputs(naomi_data, indicators, meta_area)

  val <- list(
    indicators = indicators,
    art_attendance = art_attendance,
    meta_area = meta_area,
    meta_age_group = meta_age_group,
    meta_period = meta_period,
    meta_indicator = meta_indicator,
    fit = fit,
    inputs_outputs = inputs_outputs
  )

  class(val) <- "naomi_output"

  val
}

#' Version of extract_indicators for T1 only
local_extract_indicators_naomi_simple <- function(naomi_fit, naomi_mf, na.rm = FALSE) {

  get_est <- function(varname,
                      indicator,
                      calendar_quarter,
                      mf = naomi_mf$mf_out) {
    v <- dplyr::mutate(
                  mf,
                  calendar_quarter = calendar_quarter,
                  indicator = indicator)

    tryCatch(
      if(!is.null(naomi_fit$sample)) {
        v <- naomi:::add_stats(v, naomi_fit$mode[[varname]], naomi_fit$sample[[varname]], na.rm = na.rm)
      } else {
        v <- naomi:::add_stats(v, naomi_fit$mode[[varname]], na.rm = na.rm)
      },
      "error" = function(e) {
        stop(paste0("Error simulating output for indicator:", varname, ". Please contact support for troubleshooting."))
      })

    v
  }

  indicators_t1 <- c("population_t1_out" = "population",
                     "rho_t1_out" = "prevalence",
                     "plhiv_t1_out" = "plhiv",
                     "alpha_t1_out" = "art_coverage",
                     "artnum_t1_out" = "art_current_residents",
                     "artattend_t1_out" = "art_current",
                     "untreated_plhiv_num_t1_out" = "untreated_plhiv_num",
                     "plhiv_attend_t1_out" = "plhiv_attend",
                     "untreated_plhiv_attend_t1_out" = "untreated_plhiv_attend",
                     "lambda_t1_out" = "incidence",
                     "infections_t1_out" = "infections")

  indicator_est_t1 <- Map(get_est, names(indicators_t1), indicators_t1, naomi_mf$calendar_quarter1)
  indicator_est_t1 <- dplyr::bind_rows(indicator_est_t1)

  indicators_anc_t1 <- c("anc_clients_t1_out" = "anc_clients",
                         "anc_plhiv_t1_out" = "anc_plhiv",
                         "anc_already_art_t1_out" = "anc_already_art",
                         "anc_art_new_t1_out" = "anc_art_new",
                         "anc_known_pos_t1_out" = "anc_known_pos",
                         "anc_tested_pos_t1_out" = "anc_tested_pos",
                         "anc_tested_neg_t1_out" = "anc_tested_neg",
                         "anc_rho_t1_out" = "anc_prevalence",
                         "anc_alpha_t1_out" = "anc_art_coverage")

  indicator_anc_est_t1 <- Map(get_est, names(indicators_anc_t1), indicators_anc_t1,
                              naomi_mf$calendar_quarter1, list(naomi_mf$mf_anc_out))
  indicator_anc_est_t1 <- dplyr::bind_rows(indicator_anc_est_t1)

  mf_anc_out <- naomi_mf$mf_areas %>%
    dplyr::transmute(area_id,
                     sex = "female",
                     age_group = "Y015_049")

  out <- dplyr::bind_rows(
                  indicator_est_t1,
                  indicator_anc_est_t1
  )

  dplyr::select(out, names(naomi_mf$mf_out),
                calendar_quarter, indicator, mean, se, median, mode, lower, upper)
}

#' Version of extract_art_attendance for T1 only
local_extract_art_attendance_naomi_simple <- function(naomi_fit, naomi_mf, na.rm = FALSE) {

  mode <- naomi_fit$mode

  mfout <- naomi_mf$mf_out %>%
    dplyr::mutate(out_idx = dplyr::row_number())

  v <- naomi_mf$mf_artattend %>%
    dplyr::select(reside_area_id, attend_area_id) %>%
    dplyr::mutate(sex = "both",
                  age_group = "Y000_999") %>%
    dplyr::left_join(
             dplyr::rename(mfout, reside_area_id = area_id, reside_out_idx = out_idx),
             by = c("reside_area_id", "sex", "age_group")
           ) %>%
    dplyr::left_join(
             dplyr::rename(mfout, attend_area_id = area_id, attend_out_idx = out_idx),
             by = c("attend_area_id", "sex", "age_group")
           )

  if(!is.null(mode)) {
    m_artattend_ij_t1 <- mode$artattend_ij_t1_out
    m_artnum_reside_t1 <- mode$artnum_t1_out[v$reside_out_idx]
    m_artnum_attend_t1 <- mode$artattend_t1_out[v$attend_out_idx]
  } else {
    m_artattend_ij_t1 <- NULL
  }

  if(!is.null(m_artattend_ij_t1)) {
    m_prop_residents_t1 <- m_artattend_ij_t1 / m_artnum_reside_t1
    m_prop_attendees_t1 <- m_artattend_ij_t1 / m_artnum_attend_t1
  } else {
    m_prop_residents_t1 <- NULL
    m_prop_attendees_t1 <- NULL
  }

  if(!is.null(naomi_fit$sample)) {

    s_artattend_ij_t1 <- naomi_fit$sample$artattend_ij_t1_out
    s_artnum_reside_t1 <- naomi_fit$sample$artnum_t1_out[v$reside_out_idx, ]
    s_artnum_attend_t1 <- naomi_fit$sample$artattend_t1_out[v$attend_out_idx, ]
  } else {
    s_artattend_ij_t1 <- NULL
  }

  if(!is.null(s_artattend_ij_t1)) {
    s_prop_residents_t1 <- s_artattend_ij_t1 / s_artnum_reside_t1
    s_prop_attendees_t1 <- s_artattend_ij_t1 / s_artnum_attend_t1
  } else {
    s_prop_residents_t1 <- NULL
    s_prop_attendees_t1 <- NULL
  }

  v$reside_out_idx <- NULL
  v$attend_out_idx <- NULL

  v_t1 <- dplyr::mutate(v, calendar_quarter = naomi_mf$calendar_quarter1)
  v_t1 <- naomi:::add_stats(v_t1, m_artattend_ij_t1, s_artattend_ij_t1, "artnum_", na.rm = na.rm)
  v_t1 <- naomi:::add_stats(v_t1, m_prop_residents_t1, s_prop_residents_t1, "prop_residents_", na.rm = na.rm)
  v_t1 <- naomi:::add_stats(v_t1, m_prop_attendees_t1, s_prop_attendees_t1, "prop_attendees_", na.rm = na.rm)

  dplyr::bind_rows(v_t1)
}

#' Return a quadrature rule at which to evaluate the Laplace marginal
#'
#' The grid is approximate in that it uses the Gaussian approximation to the mode
#' and standard deviation rather than the "true" mode and standard deviation which
#' one might calculate directly using the Laplace marginal. Most of the time these
#' should be very similar (Gaussian approximations tend to be good at capturing the
#' first two moments) so it should not matter.
#'
#' @param modeandhessian The row of `modesandhessians` containing the node which
#' is the mode of the Laplace approximation, or alternatively just the node which
#' has the highest log posterior evaluation.
#' @param i The index of the latent field to choose
#' @param k The number of AGHQ grid points to choose
create_approx_grid <- function(modeandhessian, i, k = 7) {
  mode <- modeandhessian[["mode"]][[1]]
  mode_i <- mode[i]
  H <- modeandhessian[["H"]][[1]]
  var_i <- diag(solve(H))[i]
  # LL <- Cholesky(H, LDL = FALSE)
  # var_i <- (colSums(solve(expand(LL)$L)^2))[i]
  gg <- mvQuad::createNIGrid(dim = 1, type = "GHe", level = k) # Create Gauss-Hermite quadrature
  mvQuad::rescale(gg, m = mode_i, C = var_i) # Adapt to mode_i and sd_i
  return(gg)
}

logSumExpWeights <- function(lp, w) {
  matrixStats::logSumExp(log(w) + lp)
}

#' For sparse grids, some of the weights are negative. This breaks the logSumExp()
#' approach unless modifications are made. The workaround is to split the sum into
#' cases when the weights are positive and cases when the weights are negative. In
#' particular, suppose that not all parts of some vector `x = (x1, ..., xm)` are
#' positive. Call the positive parts `xP` and the negative parts `xN`. Then
#' `sum(x)` is the difference between `sum(|xP|)` and `sum(|xN|)`
logSumExpNegWeights <- function(lp, w) {
  logDiffExp(
    lp1 = matrixStats::logSumExp(log(w[w > 0]) + lp[w > 0]),
    lp2 = matrixStats::logSumExp(log(-w[w < 0]) + lp[w < 0])
  )
}

#' Where we make use of calculating `log(exp(lp1) - exp(lp2))` by the following trick:
#' `log(exp(lp1) - exp(lp2)) =`
#' `log(exp(lp2) * [exp(lp1) / exp(lp2)] - exp(lp2) [exp(lp2) / exp(lp2)]) =`
#' `log(exp(lp2) * [exp(lp1 - lp2) - 1]) =`
#' `lp2 + log[exp(lp1 - lp2) - 1] =`
#' `lp2 + log(expm1(lp1 - lp2))`
logDiffExp <- function(lp1, lp2) {
  if(length(lp2) == 0) return(lp1)
  if(lp2 > lp1) {
    warning(
      "Error: the output of this function would be negative. As probabilities can't be negative, this likely isn't what you want. Returning an NA."
    )
    return(NA)
  }
  return(lp2 + log(expm1(lp1 - lp2)))
}

#' A rough unit test that the logSumExpNegWeights function does as it is intended
#' to do -- that is the logarithm of a weighted sum
stopifnot(abs(logSumExpNegWeights(lp = c(log(0.5), log(0.1)), w = c(1, -0.1)) - log(0.5 * 1 + 0.1 * -0.1)) < 10e-15)

#' Lagrange polynomial interpolant of the marginal posterior
#'
#' @param nodes Set of input values
#' @lps Log-probabilities at nodes
plot_marginal_spline <- function(nodes, lps, lower = min(nodes) - 1, upper = max(nodes) + 1) {
  ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
  interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  finegrid <- seq(lower, upper, length.out = 1000)
  df <- data.frame(x = finegrid, y = exp(interpolant(finegrid)))

  ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    theme_minimal() +
    labs(x = "x", y = "Posterior")
}

trapezoid_rule_log <- function(x, spacing) {
  w <- rep(spacing, length(x))
  w[1] <- w[1] / 2
  w[length(x)] <- w[length(x)] / 2
  matrixStats::logSumExp(log(w) + x)
}

compute_pdf_and_cdf <- function(nodes, lps, method = "auto", normalise = FALSE) {
  k <- length(nodes)
  if(k >= 4) method <- "spline"
  if(k < 4) method <- "polynomial"

  rn <- range(nodes)
  rnl <- diff(rn)
  min <- min(rn) - rnl / 2
  max <- max(rn) + rnl / 2

  if(method == "spline") {
    ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
    interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  }

  if(method == "polynomial") {
    interpolant <- as.function(polynom::poly.calc(x = nodes, y = lps))
  }

  finegrid <- seq(min, max, length.out = 1000)
  lps <- interpolant(finegrid)

  if(normalise) {
    #' Make sure that the log probabilities produce a normalised PDF
    logC <- trapezoid_rule_log(lps, spacing = finegrid[2] - finegrid[1])
    lps <- lps - logC
  }

  df <- data.frame(
    x = finegrid,
    pdf = exp(lps),
    cdf = cumsum(exp(lps)) * c(0, diff(finegrid))
  )

  return(df)
}

sample_cdf <- function(df, M) {
  q <- stats::runif(M)
  samples <- numeric(M)
  for(i in 1:M) samples[i] <- df$x[max(which(df$cdf < q[i]))]
  return(samples)
}
