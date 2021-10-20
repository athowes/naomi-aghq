tmb_summary <- function(sd_out) {
  summary(sd_out) %>%
    as.data.frame() %>%
    #' Warning! Having an issue here with rownames, needs to be fixed
    tibble::rownames_to_column() %>%
    mutate(rowname = str_replace(
      rowname,
      pattern = "\\.[0-9]+",
      replacement = paste0("[", as.character(as.numeric(str_extract(rowname, "[0-9]+")) + 1), "]"))
    ) %>%
    rename(
      "parameter" = "rowname",
      "mean" = "Estimate",
      "sd" = "Std. Error"
    ) %>%
    mutate(method = "TMB")
}

tmbstan_summary <- function(fit) {
  summary(fit)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter") %>%
    select(parameter, mean, sd) %>%
    mutate(method = "tmbstan")
}

aghq_summary <- function(quad) {
  aghq_hyper <- summary(quad)$aghqsummary$summarytable %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    rename(parameter = rowname) %>%
    select(parameter, mean, sd)

  aghq_random <- summary(quad)$randomeffectsummary %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    #' Warning: this rownumber as index only works when you have one random effect!
    mutate(variable = paste0(variable, "[", rowname, "]")) %>%
    rename(parameter = variable) %>%
    select(parameter, mean, sd)

  bind_rows(aghq_hyper, aghq_random) %>%
    mutate(method = "aghq")
}

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

  out[["comparison_results"]] <- df

  return(out)
}
