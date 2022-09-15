#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

sim_data <- readRDS("depends/sim_data.rds")
data <- sim_data[[1]] #' Only need once instance of the simulated data

compile("model1index.cpp")
dyn.load(dynlib("model1index"))

dat <- list(n = data$n, y_prev = data$y_prev, m_prev = data$m_prev)

prepare_index <- function(dat, i) {
  dat[["y_prev_i"]] <- dat$y_prev[i]
  dat[["y_prev_minus_i"]] <- dat$y_prev[-i]
  dat[["m_prev_i"]] <- dat$m_prev[i]
  dat[["m_prev_minus_i"]] <- dat$m_prev[-i]
  dat[c("n", "y_prev_i", "y_prev_minus_i", "m_prev_i", "m_prev_minus_i")]
}

#' Starting parameters for TMB are the same for every index of the loop
param <- list(
  beta_prev = 0,
  phi_prev_i = 0,
  phi_prev_minus_i = rep(0, data$n - 1),
  log_sigma_phi_prev = 0
)

#' Loop over doing the Gaussian approximation for every index i
mean <- list()
prec <- list()

for(i in 1:dat$n) {
  dat_i <- prepare_index(dat, i)

  #' random are integrated out with a Laplace approximation
  obj <- MakeADFun(
    data = dat1,
    parameters = param,
    random = "phi_prev_minus_i", # Random effects to be integrated out
    DLL = "model1index"
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

  mean[[i]] <- 1
  prec[[i]] <- sd_out$jointPrecision
}

mean
prec

pdf("prec.pdf", h = 5, w = 6.25)

lapply(prec, image)

dev.off()
