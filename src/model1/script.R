#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("model1")
# setwd("src/model1")

set.seed(1)

compile("model1.cpp")
dyn.load(dynlib("model1"))

#' Simulate data
n <- 5
m_prev <- rep(20, n)
beta_rho <- -1
tau_phi_rho <- 1
logit_rho <- beta_0 + rnorm(n, 0, 1 / sqrt(tau_phi))
rho <- plogis(logit_rho)
y_prev <- rbinom(n, m, rho)

dat <- list(n = n, y_prev = y_prev, m_prev = m_prev)

#' TMB

param <- list(
  beta_rho = 0,
  phi_rho = rep(0, n),
  log_sigma_phi_rho = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = "phi_rho", # Random effects to be integrated out
  DLL = "model1"
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
  startingvalue = c(param$beta_rho, param$log_sigma_phi_rho)
)

#' Comparison

tmb <- as.vector(t(data.frame(sd_out$par.fixed[1:2], sqrt(diag(sd_out$cov.fixed)[1:2]))))
tmbstan <- as.vector(t(summary(fit)$summary[c("beta_rho", "log_sigma_phi_rho"), c(1, 3)]))
aghq <- as.vector(t(summary(quad)$summarytable[1:2, c(1, 4)]))

df <- cbind(tmb, tmbstan, aghq) %>%
  as.data.frame() %>%
  mutate(type = gl(2, 1, 4, labels = c("Mean", "SD")))

saveRDS(df, "comparison-results.rds")
