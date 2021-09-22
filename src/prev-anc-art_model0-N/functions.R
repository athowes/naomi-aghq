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
  tmb <- as.vector(t(data.frame(sd_out$par.fixed[1], sqrt(sd_out$cov.fixed[1]))))
  tmbstan <- as.vector(t(summary(fit)$summary[1, c(1, 3)]))

  df <- cbind(tmb, tmbstan) %>%
    as.data.frame() %>%
    mutate(
      type = gl(2, 1, 2, labels = c("Mean", "SD")),
      parameter = rep(names(sd_out$par.fixed), each = 2)
    )

  out[["comparison_results"]] <- df

  return(out)
}
