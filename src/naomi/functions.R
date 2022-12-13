#' A local version of naomi::make_tmb_obj
local_make_tmb_obj <- function(data, par_init, calc_outputs = 0L, inner_verbose, progress) {
  # Begin expose naomi:::make_tmb_obj
  # https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L496
  data$calc_outputs <- as.integer(calc_outputs)

  obj <- TMB::MakeADFun(
    data = data,
    parameters = par,
    DLL = "naomi",
    silent = !inner_verbose,
    random = c("beta_rho", "beta_alpha",
               "beta_alpha_t2", "beta_lambda", "beta_asfr", "beta_anc_rho",
               "beta_anc_alpha", "beta_anc_rho_t2", "beta_anc_alpha_t2",
               "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs", "u_rho_a",
               "u_rho_as", "u_rho_xa", "u_alpha_x", "us_alpha_x",
               "u_alpha_xs", "us_alpha_xs", "u_alpha_a", "u_alpha_as",
               "u_alpha_xt", "u_alpha_xa", "u_alpha_xat", "u_alpha_xst",
               "ui_lambda_x", "logit_nu_raw", "ui_asfr_x", "ui_anc_rho_x",
               "ui_anc_alpha_x", "ui_anc_rho_xt", "ui_anc_alpha_xt",
               "log_or_gamma", "log_or_gamma_t1t2")
  )

  if (!is.null(progress)) {
    obj$fn <- report_progress(obj$fn, progress)
  }

  obj
  # End expose naomi:::make_tmb_obj
}

#' A local version of naomi::fit_tmb
local_fit_tmb <- function(tmb_input, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL) {
  # Begin expose naomi::fit_tmb
  # https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L557
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))

  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress)
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
  objout <- naomi:::make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 1L, inner_verbose)
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
  hess <- sdreport_joint_precision(fit$obj, fit$par.fixed)
  if (verbose) print("Inverting precision for joint covariance")
  cov <- solve(hess)
  if (verbose) print("Drawing sample")
  smp <- mvtnorm::rmvnorm(nsample, fit$par.full, cov)
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

#' Inference for the Naomi model using aghq
fit_aghq <- function(tmb_inputs, ...) {
  stopifnot(inherits(tmb_input, "naomi_tmb_input"))
  obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress)
  quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, ...)
  quad
}
