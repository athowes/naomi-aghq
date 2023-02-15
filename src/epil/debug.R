#' The adam code in epil is right, but the adam code in naomi-simple is wrong
#' Approach to debugging this is to copy the naomi-simple adam code into epil
#' Guess that something should go wrong, and we can debug from there

#' Run epil.Rmd up to generation of laplace_x (line 626)

#' Expose xi_laplace_marginal(i = 1, opt_theta = opt, finegrid = fine_grid)
i <- 1
opt_theta <- opt
finegrid <- fine_grid

dat$toggle <- index_controller[i, "toggle"]
dat$i <- index_controller[i, "i"]

theta_names <- unique(names(opt_theta$par))
map_fixed_theta <- list()
param_fixed_theta <- param

for(theta in theta_names) {
  map_fixed_theta[[theta]] <- rep(factor(NA), length(param[[theta]]))
  param_fixed_theta[[theta]] <- opt_theta$par[names(opt_theta$par) == theta]
}

#' We fix the beta parameter to zero, as we're using beta_i and beta_minus_i
if(dat$toggle == 0) {
  param_fixed_theta[["beta_i"]] <- param_fixed_theta$beta[i]
  param_fixed_theta[["beta_minus_i"]] <- param_fixed_theta$beta[-i]
  map_fixed_theta[["beta"]] <- rep(factor(NA), length(param[["beta"]]))
}

get_random <- function(dat) {
  if(dat$toggle == 0) return(c("beta_minus_i", "nu", "epsilon"))
  if(dat$toggle == 1) return(c("beta", "nu_minus_i", "epsilon"))
  if(dat$toggle == 2) return(c("beta", "nu", "epsilon_minus_i"))
}

obj_fixed_theta <- MakeADFun(
  data = dat,
  parameters = param_fixed_theta,
  random = get_random(dat),
  DLL = "epil_index",
  silent = TRUE,
  map = map_fixed_theta
)

#' Old way (right)
quad <- aghq::aghq(
  ff = obj_fixed_theta,
  k = 3,
  startingvalue = 0,
  control = aghq::default_control_tmb()
)

pdf_and_cdf <- aghq::compute_pdf_and_cdf(quad, finegrid = fine_grid)[[1]]

#' New way (wrong?)
spline_nodes <- function(modeandhessian, i, k = 7) {
  mode <- modeandhessian[["mode"]][[1]]
  mode_i <- mode[i]
  H <- modeandhessian[["H"]][[1]]
  var_i <- diag(solve(H))[i]
  # LL <- Cholesky(H, LDL = FALSE)
  # var_i <- (colSums(solve(expand(LL)$L)^2))[i]

  #' Create Gauss-Hermite quadrature
  gg <- mvQuad::createNIGrid(dim = 1, type = "GHe", level = k)

  #' Adapt to mode_i and sd_i
  mvQuad::rescale(gg, m = mode_i, C = var_i)

  #' Return the set of x input values
  mvQuad::getNodes(gg)
}

#' fit11 is the quad object for hypers
#' Check the weights from aghq k = 1
aghq1 <- aghq::marginal_laplace_tmb(obj, k = 1, startingvalue = c(param$l_tau_epsilon, param$l_tau_nu))
aghq1$normalized_posterior$nodesandweights

theta_mode_location <- which.max(aghq1$normalized_posterior$nodesandweights$logpost_normalized)
modeandhessian <- aghq1$modesandhessians[theta_mode_location, ]

#' These seem visually reasonable (as compared with plot from notebook) -- i.e around 1.6
(nodes <- spline_nodes(modeandhessian, 1, k = 7))

laplace_marginal <- function(x) {
  lp <- as.numeric(- obj_fixed_theta$fn(x)) + log(aghq1$normalized_posterior$nodesandweights$weights)
  return(lp - aghq1$normalized_posterior$lognormconst)
}

laplace_marginal_old <- function(x) {
  as.numeric(- obj_fixed_theta$fn(x)) - quad$normalized_posterior$lognormconst
}

lps <- purrr::map_dbl(nodes, laplace_marginal)
lps_old <- purrr::map_dbl(nodes, laplace_marginal_old)

trapezoid_rule_log <- function(x, spacing) {
  w <- rep(spacing, length(x))
  w[1] <- w[1] / 2
  w[length(x)] <- w[length(x)] / 2
  matrixStats::logSumExp(log(w) + x)
}

compute_pdf_and_cdf <- function(nodes, lps, method = "auto", normalise = TRUE) {
  k <- length(nodes)
  if(k >= 4) method <- "spline"
  if(k < 4) method <- "polynomial"

  if(method == "spline") {
    ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
    interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  }

  if(method == "polynomial") {
    interpolant <- as.function(polynom::poly.calc(x = nodes, y = lps))
  }

  lps <- interpolant(finegrid)

  if(normalise) {
    #' Make sure that the log probabilities produce a normalised PDF
    logC <- trapezoid_rule_log(lps, spacing = finegrid[2] - finegrid[1])
    lps <- lps - logC
  }

  df <- data.frame(
    theta = finegrid,
    pdf = exp(lps),
    cdf = cumsum(exp(lps)) * c(0, diff(finegrid))
  )

  return(df)
}

pdf_and_cdf_old <- compute_pdf_and_cdf(nodes, lps_old, normalise = FALSE)
pdf_and_cdf_new <- compute_pdf_and_cdf(nodes, lps, normalise = FALSE)

pdf_and_cdf_old$type <- "With lognormconst from aghq::aghq(obj_fixed_theta, ...)"
pdf_and_cdf_new$type <- "With lognromconst from aghq::aghq(obj, ...)"
pdf_and_cdf$type <- "aghq::aghq(obj_fixed_theta, ...)"

pac <- bind_rows(pdf_and_cdf_old, pdf_and_cdf_new, pdf_and_cdf)

ggplot(pac, aes(x = theta, y = cdf, col = type)) +
  geom_line(size = 1, alpha = 0.8, position = position_dodge(width = 0.01)) +
  lims(x = c(1, 2.5)) +
  labs(x = "beta[1]", y = "CDF", col = "Method") +
  theme_minimal()

#' These are pre-adaption
calc_weight <- function(d) {
  (2 * pi)^{-d/2}
}

plot(calc_weight(1:5))
