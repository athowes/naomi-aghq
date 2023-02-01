#' Start of if(adam) loop

#' Try with k = 2 and sparse grid
n_hyper <- 24
k <- 2

sparse_grid <- mvQuad::createNIGrid(n_hyper, "GHe", k, "sparse")

control <- aghq::default_control_tmb()
control$method_summaries <- "correct"
control$ndConstruction <- "sparse"

#' Start expose fit_aghq(tmb_inputs, k = k, basegrid = sparse_grid, control = control)
tmb_input <- tmb_inputs
k <- k
basegrid <- sparse_grid
control <- control
inner_verbose <- FALSE
progress <- NULL
map <- NULL
DLL <- "naomi_simple"

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL)

#' Start expose aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, ...)
ff <- obj
startingvalue <- obj$par
transformation <- aghq::default_transformation()
optresults <- NULL
basegrid <- basegrid
control <- control

validate_control(control, type = "tmb")
validate_transformation(transformation)

transformation <- make_transformation(transformation)

thetanames <- NULL

if (exists("par", ff)) thetanames <- make.unique(names(ff$par), sep = "")

if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr, theta, method = "Richardson")
}

quad <- aghq(
  ff = ff, k = k, transformation = transformation,
  startingvalue = startingvalue, optresults = optresults,
  basegrid = basegrid, control = control
)

#' Save environment here for safety / time saving
save.image(file = "adam-laplace-dev-72.rdata")
load(file = "adam-laplace-dev-72.rdata")

distinctthetas <- quad$normalized_posterior$nodesandweights[, grep("theta", colnames(quad$normalized_posterior$nodesandweights))]

if (!is.data.frame(distinctthetas)) distinctthetas <- data.frame(theta1 = distinctthetas)
modesandhessians <- distinctthetas

if (is.null(thetanames)) {
  thetanames <- colnames(distinctthetas)
} else {
  colnames(modesandhessians) <- thetanames
  colnames(quad$normalized_posterior$nodesandweights)[grep("theta", colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
}

modesandhessians$mode <- vector(mode = "list", length = nrow(distinctthetas))
modesandhessians$H <- vector(mode = "list", length = nrow(distinctthetas))

for (i in 1:nrow(distinctthetas)) {
  theta <- as.numeric(modesandhessians[i, thetanames])
  ff$fn(theta)
  mm <- ff$env$last.par
  modesandhessians[i, "mode"] <- list(list(mm[ff$env$random]))
  H <- ff$env$spHess(mm, random = TRUE)
  H <- rlang::duplicate(H)
  modesandhessians[i, "H"] <- list(list(H))
}

modesandhessians$weights <- quad$normalized_posterior$nodesandweights$weights

quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")
quad$obj <- obj

#' modesandhessians has 1 column for each hyper, and then 1 for the mode of the latents, 1 for the Hessian of the latents, and 1 for the weights
stopifnot(ncol(modesandhessians) == length(thetanames) + 3)

#' Let's begin by just trying to get the Laplace marginals working for one parameter
#' Baby steps!

compile("naomi_simple_beta_rho_index.cpp")
dyn.load(dynlib("naomi_simple_beta_rho_index"))

#' I have adapted naomi_simple_beta_rho_index to have DATA_INTEGER(i) which will
#' tell the template which index of i we are going to be Laplace approximating.

#' Create a new data object, so as to be a little less confusing
tmb_inputs_simple_i <- tmb_inputs_simple

#' Just put it to the first index for now
tmb_inputs_simple_i$data$i <- 1

#' Change the initialisation values to remove beta_rho
tmb_inputs_simple_i$par_init$beta_rho <- NULL
tmb_inputs_simple_i$par_init$beta_rho_minus_i <- 0
tmb_inputs_simple_i$par_init$beta_rho_i <- 0

#' Expose local_make_tmb_obj
data <- tmb_inputs_simple_i$data
par <- tmb_inputs_simple_i$par_init
calc_outputs <- 0L
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL
map <- NULL

#' Change the DLL to the beta_rho_index version
DLL <- "naomi_simple_beta_rho_index"

data$calc_outputs <- as.integer(calc_outputs)

#' Same as naomi_simple but with beta_rho_minus_i rather than beta_rho
integrate_out <- c(
  "beta_rho_minus_i", "beta_alpha", "beta_lambda", "beta_anc_rho",
  "beta_anc_alpha", "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs",
  "u_rho_a", "u_rho_as", "u_rho_xa", "u_alpha_x", "us_alpha_x",
  "u_alpha_xs", "us_alpha_xs", "u_alpha_a", "u_alpha_as", "u_alpha_xa",
  "ui_lambda_x", "ui_anc_rho_x", "ui_anc_alpha_x", "log_or_gamma"
)

#' Set-up so that the initial random effects are chosen to be the parameter of the
#' latest likelihood evaluation. I will later manipulate this option to put whatever
#' random effects I'd like to start with as the "latest" likelihood evaluation
obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = DLL,
  silent = !inner_verbose,
  random = integrate_out,
  map = map,
  random.start = expression(last.par[random])
)

if (!is.null(progress)) {
  obj$fn <- naomi:::report_progress(obj$fn, progress)
}

#' A function to choose the points at which to evaluate the Laplace marginal
#'
#' @param modeandhessian The row of `modesandhessians` containing the node which
#' is the mode of the Laplace approximation, or alternatively just the node which
#' has the highest log posterior evaluation.
#' @param i The index of the latent field to choose
#' @param k The number of AGHQ grid points to choose
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

theta_mode_location <- which.max(quad$normalized_posterior$nodesandweights$logpost_normalized)
modeandhessian <- modesandhessians[theta_mode_location, ]

#' The set of input values that we'd like to calculate the log-probability at
(nodes <- spline_nodes(modeandhessian, 1, k = 7))

#' Check the order of parameters in obj
obj$par

#' For sparse grids, some of the weights are negative. This breaks the logSumExp()
#' approach unless modifications are made. The workaround is to split the sum into
#' cases when the weights are positive and cases when the weights are negative. In
#' particular, suppose that not all parts of some vector `x = (x1, ..., xm)` are
#' positive. Call the positive parts `xP` and the negative parts `xN`. Then
#' `sum(x)` is the difference between `sum(|xP|)` and `sum(|xN|)`
logSumExpWeights <- function(lp, w) {
  logDiffExp(
    lp1 = matrixStats::logSumExp(log(w[w > 0]) + lp[w > 0]),
    lp2 = matrixStats::logSumExp(log(-w[w < 0]) + lp[w < 0])
  )
}

#' Where we make use of calculating `log(exp(lp1) - exp(lp2))` by the following trick:
#' `log(exp(lp1) - exp(lp2)) =`
#' `log(exp(lp2) * [exp(lp1) / exp(lp2)] - exp(lp2) [exp(lp2) / exp(lp2)]) =`
#' `log(exp(lp2) * [exp(lp1 - lp2) - 1]) =`
#' `lp2 + log[exp(lp1 - lp2) - 1] =`
#' `lp2 + log(expm1(lp1 - lp2))`
logDiffExp <- function(lp1, lp2) {
  if(length(lp2) == 0) return(lp1)
  if(lp2 > lp1) {
    warning(
      "Error: the output of this function would be negative. As probabilities can't be negative, this likely isn't what you want. Returning an NA."
    )
    return(NA)
  }
  return(lp2 + log(expm1(lp1 - lp2)))
}

#' A rough unit test that the logSumExpWeights function does as it is intended
#' to do -- that is the logarithm of a weighted sum
stopifnot(abs(logSumExpWeights(lp = c(log(0.5), log(0.1)), w = c(1, -0.1)) - log(0.5 * 1 + 0.1 * -0.1)) < 10e-15)

laplace_marginal <- function(x) {
  lp <- vector(mode = "numeric", length = nrow(modesandhessians))
  random <- obj$env$random   #' Indices of random effects

  for(z in 1:nrow(modesandhessians)) {
    theta <- as.numeric(modesandhessians[z, thetanames])
    #' Get the mode of the random effects at the hyperparameter node location
    #' and set things such that obj$fn is initialised there using last.par.
    #' The -1 removes the element of beta_rho that we are Laplace approximating
    #' and has been moved out of the random into the hyper
    mode <- modesandhessians[z, "mode"][[1]][-1]
    obj$env$last.par[random] <- mode
    lp[z] <- as.numeric(- obj$fn(c(x, theta)))
  }

  logSumExpWeights(lp, w = modesandhessians$weights) - quad$normalized_posterior$lognormconst
}

#' The log-probabilities at the set of input values. I'm also timing how long it
#' takes, as this will determine how viable this method is when moving to all
#' 500 or so marginals
lps <- vector(mode = "numeric", length = length(nodes))
starts <- vector(mode = "numeric", length = length(nodes))
ends <- vector(mode = "numeric", length = length(nodes))
times <- vector(mode = "numeric", length = length(nodes))

for(i in seq_along(nodes)) {
  starts[i] <- Sys.time()
  lps[i] <- laplace_marginal(nodes[i])
  ends[i] <- Sys.time()
  times[i] <- ends[i] - starts[i]
}

#' Lagrange polynomial interpolant of the marginal posterior
#'
#' @param nodes Set of input values
#' @lps Log-probabilities at nodes
plot_marginal_spline <- function(nodes, lps) {
  ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
  interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  finegrid <- seq(-5, 5, by = 0.1)
  df <- data.frame(x = finegrid, y = exp(interpolant(finegrid)))

  ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    theme_minimal() +
    labs(x = "x", y = "Posterior")
}

plot_marginal_spline(nodes, lps)

#' Comparison to inference results from other methods
H <- modeandhessian[["H"]][[1]]
var_i <- diag(solve(H))[i]

beta_rho <- readRDS("depends/beta_rho.rds")

#' This is what the posterior marginal should look like (from aghq, TMB, tmbstan)
plot <- beta_rho %>%
  ggplot(aes(x = samples, fill = method, col = method)) +
  geom_density(aes(y = after_stat(density)), alpha = 0.05, position = "identity") +
  theme_minimal() +
  labs(x = "beta_rho[1]", y = "Density", col = "Method") +
  scale_color_manual(values = multi.utils::cbpalette()) +
  scale_fill_manual(values = multi.utils::cbpalette()) +
  guides(fill = FALSE)

plot +
  geom_point(
    data = data.frame(x = nodes, y = exp(lps)) %>%
      filter(!is.na(lps)),
    aes(x = x, y = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    x = -10, y = 0.35,
    label = "Black points represent Laplace marginal evaluations",
    inherit.aes = FALSE
  )

#' Trying to understand the node which give NA results
laplace_marginal(nodes[3])

x <- nodes[3]

lp <- vector(mode = "numeric", length = nrow(modesandhessians))
random <- obj$env$random   #' Indices of random effects

for(z in 1:nrow(modesandhessians)) {
  theta <- as.numeric(modesandhessians[z, thetanames])
  mode <- modesandhessians[z, "mode"][[1]][-1]
  obj$env$last.par[random] <- mode
  lp[z] <- as.numeric(- obj$fn(c(x, theta)))
}

lp_normalised <- lp - quad$normalized_posterior$lognormconst

logSumExpWeights(lp_normalised, w = modesandhessians$weights)

lp1 <- matrixStats::logSumExp(log(modesandhessians$weights[modesandhessians$weights > 0]) + lp_normalised[modesandhessians$weights > 0])
lp2 <- matrixStats::logSumExp(log(-modesandhessians$weights[modesandhessians$weights < 0]) + lp_normalised[modesandhessians$weights < 0])

data.frame(x = modesandhessians$weights, y = lp_normalised) %>%
  ggplot(aes(x = x, y = y)) +
    geom_point() +
    labs(x = "Weight", y = "Log-posterior") +
    theme_minimal()

#' Aiming to get the marginal standard deviation without inverting the whole Hessian
H <- modeandhessian[["H"]][[1]]
LL <- Cholesky(H, LDL = FALSE)
dd1 <- diag(solve(H))
dd2 <- colSums(solve(expand(LL)$L)^2)
sum(abs(dd1 - dd2)) #' Different answers?
