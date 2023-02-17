#' Compare approaches to calculating one particular Laplace posterior marginal
#' to see if there are any differences. We choose the beta_anc_rho marginal,
#' because previously it has given problems where the PDF integrates to > 1.
#'
#' We will build on the code in k1-adam-laplace-one-x.R

quad <- readRDS("depends/out.rds")$quad

modesandhessians <- quad$modesandhessians
theta_names <- make.unique(names(quad$obj$par), sep = "")
modesandhessians$weights <- quad$normalized_posterior$nodesandweights$weights

random <- quad$obj$env$random
xnames <- names(quad$obj$env$par[random])

tmb_inputs_simple_i <- tmb_inputs_simple
x_lengths <- lengths(tmb_inputs_simple_i$par_init[unique(xnames)])
tmb_inputs_simple_i$par_init$x_minus_i <- rep(0, sum(x_lengths) - 1)
tmb_inputs_simple_i$par_init$x_i <- 0
tmb_inputs_simple_i$par_init[unique(xnames)] <- NULL

x_starts <- cumsum(x_lengths) - x_lengths
tmb_inputs_simple_i$data$x_lengths <- as.numeric(x_lengths)
tmb_inputs_simple_i$data$x_starts <- as.numeric(x_starts)
tmb_inputs_simple_i$data$i <- 7 #' "beta_anc_rho"

DLL <- "naomi_simple_x_index"
compile(paste0(DLL, ".cpp"))
dyn.load(dynlib(DLL))

data <- tmb_inputs_simple_i$data
par <- tmb_inputs_simple_i$par_init
calc_outputs <- 0L
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL
map <- NULL

data$calc_outputs <- as.integer(calc_outputs)

obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = DLL,
  silent = !inner_verbose,
  random = "x_minus_i",
  map = map,
  random.start = expression(last.par[random])
)

(nodes <- spline_nodes(modeandhessian = modesandhessians, tmb_inputs_simple_i$data$i, k = 7))

laplace_marginal <- function(x, i) {
  random <- obj$env$random
  theta <- as.numeric(modesandhessians[theta_names])
  mode <- modesandhessians[["mode"]][[1]][-i]
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
  lps[i] <- laplace_marginal(nodes[i], i = tmb_inputs_simple_i$data$i)
  ends[i] <- Sys.time()
  times[i] <- ends[i] - starts[i]
}

pdf_and_cdf <- compute_pdf_and_cdf(nodes, lps, normalise = FALSE)

#' The PDF obtained
ggplot(pdf_and_cdf, aes(x = x, y = pdf)) +
  geom_line() +
  theme_minimal()

max(pdf_and_cdf$cdf)

#' Alternative method using aghq::aghq

param_fixed_theta <- par
map_fixed_theta <- list()

for(theta in theta_names) {
  map_fixed_theta[[theta]] <- rep(factor(NA), length(par[[theta]]))
  param_fixed_theta[[theta]] <- as.numeric(quad$optresults$mode[names(quad$optresults$mode) == theta])
}

obj_fixed_theta <- TMB::MakeADFun(
  data = data,
  parameters = param_fixed_theta,
  DLL = DLL,
  silent = !inner_verbose,
  random = "x_minus_i",
  map = map_fixed_theta,
  random.start = expression(last.par[random])
)

maginal_quad <- aghq::aghq(
  ff = obj_fixed_theta,
  k = 7,
  startingvalue = 0,
  control = aghq::default_control_tmb()
)

rn <- range(nodes)
rnl <- diff(rn)
min <- min(rn) - rnl / 2
max <- max(rn) + rnl / 2
fine_grid <- seq(min, max, length.out = 1000)

pdf_and_cdf_alt <- aghq::compute_pdf_and_cdf(maginal_quad, finegrid = fine_grid)[[1]]

max(pdf_and_cdf_alt$cdf)
