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
  tmb <- as.vector(t(data.frame(sd_out$par.fixed[1:2], sqrt(diag(sd_out$cov.fixed)[1:2]))))
  tmbstan <- as.vector(t(summary(fit)$summary[c("beta_prev", "log_sigma_phi_prev"), c(1, 3)]))
  aghq <- as.vector(t(summary(quad)$summarytable[1:2, c(1, 4)]))

  df <- cbind(tmb, tmbstan, aghq) %>%
    as.data.frame() %>%
    mutate(type = gl(2, 1, 4, labels = c("Mean", "SD")))

  out[["comparison_results"]] <- df

  return(out)
}
