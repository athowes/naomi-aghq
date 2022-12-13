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
