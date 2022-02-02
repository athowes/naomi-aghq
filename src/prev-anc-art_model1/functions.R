run_model1 <- function(data) {

  dat <- list(n = data$n, y_prev = data$y_prev, m_prev = data$m_prev)

  #' TMB
  param <- list(
    beta_prev = 0,
    phi_prev = rep(0, data$n),
    log_sigma_phi_prev = 0
  )

  #' random are integrated out with a Laplace approximation
  obj <- MakeADFun(
    data = dat,
    parameters = param,
    random = "phi_prev", # Random effects to be integrated out
    DLL = "model1"
  )

  its <- 1000
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = its, trace = 0)
  )

  sd_out <- sdreport(
    obj,
    par.fixed = opt$par,
    getJointPrecision = TRUE
  )

  #' tmbstan
  fit <- tmbstan(obj = obj, chains = 4)

  out <- list(
    "mcmc_traceplots" = rstan::traceplot(fit, pars = names(obj$par), inc_warmup = TRUE),
    "mcmc_monitor" =  rstan::monitor(fit)
  )

  #' aghq
  quad <- aghq::marginal_laplace_tmb(
    obj,
    k = 3,
    startingvalue = c(param$beta_prev, param$log_sigma_phi_prev)
  )

  #' Comparison
  tmb <- tmb_summary(sd_out)
  tmbstan <- tmbstan_summary(fit)
  aghq <- aghq_summary(quad)

  df <- bind_rows(tmb, tmbstan, aghq)

  out[["comparison_results"]] <- df

  #' Kolmogorov-Smirnov test

  samples_aghq <- aghq::sample_marginal(quad, M = 1000)$samps %>%
    t() %>%
    as.data.frame() %>%
    inf.utils::replace_duplicate_colnames()

  samples_tmbstan <- as.data.frame(fit)

  samples_tmb <- sample_tmb(sd_out, obj, M = 1000) %>%
    as.data.frame() %>%
    inf.utils::replace_duplicate_colnames()

  ks_aghq_tmbstan <- lapply(colnames(samples_aghq), function(col) {
    c("parameter" = col, "ks" = get_ks(samples_aghq[, col], samples_tmbstan[, col]))
  }) %>%
    bind_rows() %>%
    mutate(method1 = "aghq", method2 = "tmbstan")

  ks_tmb_tmbstan <- lapply(colnames(samples_tmb), function(col) {
    c("parameter" = col, "ks" = get_ks(samples_tmb[, col], samples_tmbstan[, col]))
  }) %>%
    bind_rows() %>%
    mutate(method1 = "TMB", method2 = "tmbstan")

  out[["ks_test"]] <- bind_rows(ks_aghq_tmbstan, ks_tmb_tmbstan)

  return(out)
}
