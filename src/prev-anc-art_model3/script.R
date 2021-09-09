#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model3")
# setwd("src/prev-anc-art_model3")

set.seed(1)

compile("model3.cpp")
dyn.load(dynlib("model3"))

#' Simulate data
n <- 5
m_prev <- rep(20, n)
beta_prev <- -1
tau_phi_prev <- 1
eta_prev <- beta_prev + rnorm(n, 0, 1 / sqrt(tau_phi_prev))
rho_prev <- plogis(eta_prev)
y_prev <- rbinom(n, m_prev, rho_prev)
N_art <- rep(500, n)
beta_art <- 0
tau_phi_art <- 1
eta_art <- beta_art + rnorm(n, 0, 1 / sqrt(tau_phi_art))
alpha_art <- plogis(eta_art)
A_art <- rbinom(n, floor(N_art * rho_prev), alpha_art)

dat <- list(
  n = n,
  y_prev = y_prev,
  m_prev = m_prev,
  A_art = A_art,
  N_art = N_art
)

#' TMB

param <- list(
  beta_prev = 0,
  phi_prev = rep(0, n),
  log_sigma_phi_prev = 0,
  beta_art = 0,
  phi_art = rep(0, n),
  log_sigma_phi_art = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = c("phi_prev", "phi_art"),
  DLL = "model3"
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
      param$beta_art,
      param$log_sigma_phi_art
    )
)

#' Comparison

tmb <- as.vector(t(data.frame(sd_out$par.fixed[1:4], sqrt(diag(sd_out$cov.fixed)[1:4]))))
tmbstan <- as.vector(t(summary(fit)$summary[c("beta_prev", "log_sigma_phi_prev", "beta_art", "log_sigma_phi_art"), c(1, 3)]))
aghq <- as.vector(t(summary(quad)$summarytable[1:4, c(1, 4)]))

df <- cbind(tmb, tmbstan, aghq) %>%
  as.data.frame() %>%
  mutate(type = gl(2, 1, 8, labels = c("Mean", "SD")))

saveRDS(df, "comparison-results.rds")
