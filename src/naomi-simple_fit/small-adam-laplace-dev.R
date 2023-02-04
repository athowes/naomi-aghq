#' This development script attempts to get AGHQ with Laplace marignals working
#' for AGHQ k = 1. It's small in the sense of being simple enough that it works,
#' which I have now confirmed. The next steps here are to get this into proper
#' function.R code.

#' Going to do k = 1
#' Assume that the if(aghq) loop has been run
#' Can just get this from a saved version
quad <- readRDS("depends/out.rds")$quad

modesandhessians <- quad$modesandhessians
thetanames <- make.unique(names(quad$obj$par), sep = "")
modesandhessians$weights <- quad$normalized_posterior$nodesandweights$weights

#' Let's begin by just trying to get the Laplace marginals working for one parameter

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

#' The set of input values that we'd like to calculate the log-probability at
#' For k = 1 there is only one row in modesandhessians so we can just pass that
(nodes <- spline_nodes(modeandhessian = modesandhessians, 1, k = 7))

#' Can be quite simple because there are no weights and it's just EB
laplace_marginal <- function(x) {
  random <- obj$env$random
  theta <- as.numeric(modesandhessians[thetanames])
  mode <- modesandhessians[["mode"]][[1]][-1]
  obj$env$last.par[random] <- mode
  lp <- as.numeric(- obj$fn(c(x, theta))) + log(modesandhessians$weights)
  lp - quad$normalized_posterior$lognormconst
}

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

plot_marginal_spline(nodes, lps, lower = min(nodes) - 1, upper = max(nodes) + 1)

#' Checking it looks the same as inference from other methods
beta_rho <- readRDS("depends/beta_rho.rds")

#' This is what the posterior marginal should look like (from aghq, TMB, tmbstan)
plot <- beta_rho %>%
  ggplot(aes(x = samples, fill = method, col = method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.05, position = "identity", bins = 30) +
  theme_minimal() +
  labs(x = "beta_rho[1]", y = "Density", col = "Method") +
  scale_color_manual(values = multi.utils::cbpalette()) +
  scale_fill_manual(values = multi.utils::cbpalette()) +
  guides(fill = FALSE)

lower <- min(nodes) - 1
upper <- max(nodes) + 1
ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
finegrid <- seq(lower, upper, by = 0.1)
df <- data.frame(x = finegrid, y = exp(interpolant(finegrid)))

plot +
  geom_line(
    data = df,
    aes(x = x, y = y),
    inherit.aes = FALSE
  )
