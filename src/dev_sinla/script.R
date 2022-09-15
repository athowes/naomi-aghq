#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

sim_data <- readRDS("depends/sim_data.rds")
data <- sim_data[[1]] #' Only need once instance of the simulated data

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

#' aghq
quad <- aghq::marginal_laplace_tmb(
  obj,
  k = 3,
  startingvalue = c(param$beta_prev, param$log_sigma_phi_prev)
)

#' Comparison
tmb <- tmb_summary(sd_out)
aghq <- aghq_summary(quad)
