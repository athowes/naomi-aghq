#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("dev_sinla")
# setwd("src/dev_sinla")

sim_data <- readRDS("depends/sim_data.rds")
data <- sim_data[[1]] #' Only need once instance of the simulated data

compile("model1index.cpp")
dyn.load(dynlib("model1index"))

dat <- list(n = data$n, y_prev = data$y_prev, m_prev = data$m_prev)

prepare_dat <- function(dat, i) {
  dat[["y_prev_i"]] <- dat$y_prev[i]
  dat[["y_prev_minus_i"]] <- dat$y_prev[-i]
  dat[["m_prev_i"]] <- dat$m_prev[i]
  dat[["m_prev_minus_i"]] <- dat$m_prev[-i]
  dat[c("n", "y_prev_i", "y_prev_minus_i", "m_prev_i", "m_prev_minus_i")]
}

#' Loop over doing the Gaussian approximation for every index i
template <- list() #' To store obj

for(i in 1:dat$n) {
  dat_i <- prepare_index(dat, i)

  #' Starting parameters for TMB are the same for every index of the loop
  param_i <- list(
    beta_prev = 0,
    phi_prev_i = 0,
    phi_prev_minus_i = rep(0, data$n - 1),
    log_sigma_phi_prev = 0
  )

  #' random are integrated out with a Laplace approximation
  obj <- MakeADFun(
    data = dat_i,
    parameters = param_i,
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

  template[[i]] <- obj
}

# Create objective function (no Laplace approximation)
param <- list(
  beta_prev = 0,
  phi_prev = rep(0, data$n),
  log_sigma_phi_prev = 0
)

f <- MakeADFun(
  data = dat,
  parameters = param,
  DLL = "model1"
)

#' Evaluate objective function at param
f$fn(unlist(param))

transform_param <- function(param) {
  list(
    beta_prev = param$beta_prev,
    phi_prev_i = param$phi_prev[i],
    phi_prev_minus_i = param$phi_prev[-i],
    log_sigma_phi_prev = param$log_sigma_phi_prev
  )
}

param %>%
  transform_param() %>%
  unlist() %>%
  template[[1]]$fn()

mean <- with(template[[1]]$env, last.par[random])
Q <- template[[1]]$env$spHess(with(template[[1]]$env, last.par), random = TRUE)

#' Numerator of the LA
#' TODO

#' Denominator of the LA
sqrt(det(Q))
