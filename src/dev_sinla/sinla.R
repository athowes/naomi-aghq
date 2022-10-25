#' Clean script version of sinla.Rmd notebook
orderly::orderly_develop_start("dev_sinla")

cbpalette <- multi.utils::cbpalette()

sim_data <- readRDS("depends/sim_data.rds")
data <- sim_data[[1]]
dat <- list(n = data$n, y_prev = data$y_prev, m_prev = data$m_prev)

compile("model1.cpp")
dyn.load(dynlib("model1"))

param <- list(
  beta_prev = -2,
  phi_prev = rep(0, data$n),
  log_sigma_phi_prev = -1
)

h <- MakeADFun(
  data = dat,
  parameters = param,
  random = "phi_prev",
  DLL = "model1",
  silent = TRUE
)

#' 1: MCMC
mcmc <- tmbstan::tmbstan(h, chains = 4, iter = 5000, refresh = 0)

samples_tmbstan <- extract(mcmc, pars = "phi_prev") %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "index",
    names_prefix = "phi_prev.",
    names_transform = as.integer
  )

#' 2: Gaussian
its <- 1000

opt_theta <- nlminb(
  start = h$par,
  objective = h$fn,
  gradient = h$gr,
  control = list(iter.max = its, trace = 0)
)

mm <- h$env$last.par
mean <- mm[h$env$random]
Q <- h$env$spHess(mm, random = TRUE)
Sigma <- solve(Q)

gaussian_df <- lapply(1:length(mean), function(i) {
  x <- seq(-5, 5, length.out = 1000)
  mean_i <- as.numeric(mean[i])
  sd_i <- as.numeric(sqrt(Sigma[i, i]))
  data.frame(index = i, x = x, pdf = dnorm(x, mean = mean_i, sd = sd_i))
}) %>%
  bind_rows()

#' 3: Laplace

prepare_dat <- function(dat, i) {
  dat[["y_prev_i"]] <- dat$y_prev[i]
  dat[["y_prev_minus_i"]] <- dat$y_prev[-i]
  dat[["m_prev_i"]] <- dat$m_prev[i]
  dat[["m_prev_minus_i"]] <- dat$m_prev[-i]
  dat[c("n", "y_prev_i", "y_prev_minus_i", "m_prev_i", "m_prev_minus_i")]
}

compile("model1index.cpp")
dyn.load(dynlib("model1index"))

xi_laplace_marginal <- function(i, opt_theta) {
  dat_i <- prepare_dat(dat, i)
  map_fixed_theta <- list(beta_prev = factor(NA), log_sigma_phi_prev = factor(NA))

  param_i_fixed_theta <- list(
    beta_prev = 0,
    phi_prev_i = 0,
    phi_prev_minus_i = rep(0, data$n - 1),
    log_sigma_phi_prev = 0
  )

  for(theta in names(opt_theta$par)) param_i_fixed_theta[[theta]] <- opt_theta$par[[theta]]

  obj_fixed_theta <- MakeADFun(
    data = dat_i,
    parameters = param_i_fixed_theta,
    random = "phi_prev_minus_i",
    DLL = "model1index",
    silent = TRUE,
    map = map_fixed_theta
  )

  quad <- aghq::aghq(
    ff = obj_fixed_theta,
    k = 3,
    startingvalue = 0,
    control = aghq::default_control_tmb()
  )

  pdf_and_cdf <- aghq::compute_pdf_and_cdf(quad)[[1]]

  return(pdf_and_cdf)
}

laplace_df <- lapply(1:dat$n, xi_laplace_marginal, opt_theta = opt_theta) %>%
  bind_rows(.id = "index") %>%
  mutate(index = as.numeric(index))

#' Outputs

pdf("densities.pdf", h = 8, w = 6.25)

ggplot(samples_tmbstan, aes(x = value)) +
  geom_histogram(aes(y = ..density..), alpha = 0.8, fill = cbpalette[7]) +
  geom_line(data = gaussian_df, aes(x = x, y = pdf), col = cbpalette[1], alpha = 0.8) +
  geom_line(data = laplace_df, aes(x = theta, y = pdf), col = cbpalette[2], alpha = 0.8) +
  facet_wrap(~index) +
  theme_minimal() +
  labs(x = "phi_prev", y = "Posterior PDF")

dev.off()
