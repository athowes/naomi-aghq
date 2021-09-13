#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model1")
# setwd("src/prev-anc-art_model1")

set.seed(1)

load("depends/sim_data.rdata")

dat <- list(n = n, y_prev = y_prev, m_prev = m_prev)

#' TMB

compile("model1.cpp")
dyn.load(dynlib("model1"))

param <- list(
  beta_prev = 0,
  phi_prev = rep(0, n),
  log_sigma_phi_prev = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = "phi_prev", # Random effects to be integrated out
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

pdf("mcmc-traceplots.pdf", h = length(names(obj$par)) * 3.5, w = 8.5)

rstan::traceplot(fit, pars = names(obj$par), inc_warmup = TRUE)

dev.off()

sink(file = "mcmc-monitor.txt")

mon <- rstan::monitor(fit)
mon

sink(file = NULL)

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

saveRDS(df, "comparison-results.rds")
