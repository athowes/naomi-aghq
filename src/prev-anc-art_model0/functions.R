run_model0 <- function(data) {

  dat <- list(n = data$n, y_prev = data$y_prev, m_prev = data$m_prev)

  #' TMB
  param <- list(beta_prev = 0)

  #' random are integrated out with a Laplace approximation
  obj <- MakeADFun(
    data = dat,
    parameters = param,
    DLL = "model0"
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

  #' Comparison
  tmb <- tmb_summary(sd_out)
  tmbstan <- tmbstan_summary(fit)

  df <- bind_rows(tmb, tmbstan)

  #' Kolmogorov-Smirnov test

  samples_tmbstan <- as.data.frame(fit)

  samples_tmb <- sample_tmb(sd_out, obj, M = 1000) %>%
    as.data.frame() %>%
    inf.utils::replace_duplicate_colnames()

  ks_tmb_tmbstan <- lapply(colnames(samples_tmb), function(col) {
    ks_result <- inf.utils::ks_test(samples_tmb[, col], samples_tmbstan[, col])
    c("parameter" = col, "D" = ks_result$D, "l" = ks_result$l)
  }) %>%
    bind_rows() %>%
    mutate(method1 = "TMB", method2 = "tmbstan")

  out[["comparison_results"]] <- ks_tmb_tmbstan

  return(out)
}
