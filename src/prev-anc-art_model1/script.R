#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model1")
# setwd("src/prev-anc-art_model1")

sim_data <- readRDS("depends/sim_data.rds")

#' TMB preparation
compile("model1.cpp")
dyn.load(dynlib("model1"))

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
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  tibble::rownames_to_column("rowname") %>%
  mutate(variable = paste0(str_split(variable, pattern = "[.]")[[1]][1], "[", rowname, "]")) %>%
  select(-rowname) %>%
  tibble::column_to_rownames("variable") %>%
  t()

samples_tmbstan <- as.data.frame(fit)

#' Development starts here
sample_tmb <- function(sd_out, obj) {
  prec <- sd_out$jointPrecision
  cov <- solve(prec)
  samples <- mvtnorm::rmvnorm(n = 1000, obj$env$last.par, as.matrix(cov))
}

samples_tmb <- sample_tmb(sd_out, obj)
#' Development ends here

ks_aghq_tmbstan <- lapply(colnames(samples_aghq), function(col) {
  c("parameter" = col, "ks" = get_ks(samples_aghq[, col], samples_tmbstan[, col]))
}) %>%
  bind_rows() %>%
  mutate(method1 = "aghq", method2 = "tmbstan")

#' TODO: KS test for TMB and tmbstan
#' Requries function to sample from generic TMB model!
#' Will then bind_rows() these together before outputting below

out[["ks_test"]] <- ks_aghq_tmbstan


#' Run for each simulated dataset (using head for now for speed)
results <- map(head(sim_data), run_model1)

#' Save to artefact
saveRDS(results, "results.rds")
