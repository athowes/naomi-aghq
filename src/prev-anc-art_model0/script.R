#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("prev-anc-art_model0")
# setwd("src/prev-anc-art_model0")

set.seed(1)

load("depends/sim_data.rdata")

dat <- list(n = n, y_prev = y_prev, m_prev = m_prev)

#' TMB

compile("model0.cpp")
dyn.load(dynlib("model0"))

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

pdf("mcmc-traceplots.pdf", h = length(names(obj$par)) * 3.5, w = 8.5)

rstan::traceplot(fit, pars = names(obj$par), inc_warmup = TRUE)

dev.off()

sink(file = "mcmc-monitor.txt")

mon <- rstan::monitor(fit)
mon

sink(file = NULL)

#' aghq

#' Error in ff$env$spHess(mm, random = TRUE) : attempt to apply non-function
#' Is this because not using any random?
# quad <- aghq::marginal_laplace_tmb(
#   obj,
#   k = 3,
#   startingvalue = c(param$beta_prev)
# )

#' Comparison

tmb <- as.vector(t(data.frame(sd_out$par.fixed[1], sqrt(sd_out$cov.fixed[1]))))
tmbstan <- as.vector(t(summary(fit)$summary[1, c(1, 3)]))

df <- cbind(tmb, tmbstan) %>%
  as.data.frame() %>%
  mutate(type = gl(2, 1, 2, labels = c("Mean", "SD")))

saveRDS(df, "comparison-results.rds")
