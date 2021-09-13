#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model4")
# setwd("src/prev-anc-art_model4")

set.seed(1)

load("depends/sim_data.rdata")

dat <- list(
  n = n,
  y_prev = y_prev,
  m_prev = m_prev,
  y_anc = y_anc,
  m_anc = m_anc,
  A_art = A_art,
  N_art = N_art
)

#' TMB

compile("model4.cpp")
dyn.load(dynlib("model4"))

param <- list(
  beta_prev = 0,
  phi_prev = rep(0, n),
  log_sigma_phi_prev = 0,
  beta_anc = 0,
  b_anc = rep(0, n),
  log_sigma_b_anc = 0,
  beta_art = 0,
  phi_art = rep(0, n),
  log_sigma_phi_art = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = c("phi_prev", "b_anc", "phi_art"),
  DLL = "model4"
)

its <- 1000 # May converge before this
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

#' aghq

quad <- aghq::marginal_laplace_tmb(
  obj,
  k = 3,
  startingvalue =
    c(
      param$beta_prev,
      param$log_sigma_phi_prev,
      param$beta_anc,
      param$log_sigma_b_anc,
      param$beta_art,
      param$log_sigma_phi_art
    )
)

#' Comparison

tmb <- as.vector(t(data.frame(sd_out$par.fixed[1:6], sqrt(diag(sd_out$cov.fixed)[1:6]))))
tmbstan <- as.vector(t(summary(fit)$summary[c("beta_prev", "log_sigma_phi_prev", "beta_anc", "log_sigma_b_anc", "beta_art", "log_sigma_phi_art"), c(1, 3)]))
aghq <- as.vector(t(summary(quad)$summarytable[1:6, c(1, 6)]))

df <- cbind(tmb, tmbstan, aghq) %>%
  as.data.frame() %>%
  mutate(type = gl(2, 1, 12, labels = c("Mean", "SD")))

saveRDS(df, "comparison-results.rds")
