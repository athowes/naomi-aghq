#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("model2")
# setwd("src/model2")

set.seed(1)

compile("model2.cpp")
dyn.load(dynlib("model2"))

#' Simulate data
n <- 5
m_prev <- rep(20, n)
beta_prev <- -1
tau_phi_prev <- 1
logit_rho_prev <- beta_prev + rnorm(n, 0, 1 / sqrt(tau_phi_prev))
rho_prev <- plogis(logit_rho_prev)
y_prev <- rbinom(n, m_prev, rho_prev)
m_anc <- rep(50, n)
beta_anc <- 0.5
tau_b_anc <- 4
b_anc <- beta_anc + rnorm(n, 0,  1 / sqrt(tau_b_anc))
logit_rho_anc <- logit_rho_prev + b_anc
rho_anc <- plogis(logit_rho_anc)
y_anc <- rbinom(n, m_anc, rho_anc)

dat <- list(
  n = n,
  y_prev = y_prev,
  m_prev = m_prev,
  y_anc = y_anc,
  m_anc = m_anc
)

#' TMB

param <- list(
  beta_prev = 0,
  phi_prev = rep(0, n),
  log_sigma_phi_prev = 0,
  beta_anc = 0,
  b_anc = rep(0, n),
  log_sigma_b_anc = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = c("phi_prev", "b_anc"),
  DLL = "model2"
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
      param$log_sigma_b_anc
    )
)

#' Comparison

sd_out$par.fixed[1:4]

tmb <- as.vector(t(data.frame(sd_out$par.fixed[1:4], sqrt(diag(sd_out$cov.fixed)[1:4]))))
tmbstan <- as.vector(t(summary(fit)$summary[c("beta_prev", "log_sigma_phi_prev", "beta_anc", "log_sigma_b_anc"), c(1, 3)]))
aghq <- as.vector(t(summary(quad)$summarytable[1:4, c(1, 4)]))

df <- cbind(tmb, tmbstan, aghq) %>%
  as.data.frame() %>%
  mutate(type = gl(2, 1, 8, labels = c("Mean", "SD")))

saveRDS(df, "comparison-results.rds")
