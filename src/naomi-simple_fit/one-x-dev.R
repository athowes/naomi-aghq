#' Work on concatenating all the latent field into one x

quad <- readRDS("depends/out.rds")$quad

modesandhessians <- quad$modesandhessians
thetanames <- make.unique(names(quad$obj$par), sep = "")
modesandhessians$weights <- quad$normalized_posterior$nodesandweights$weights

random <- quad$obj$env$random
xnames <- names(quad$obj$env$par[random])

#' Replace the $par_init of tmb_inputs_simple with a concatendated version
#' which has all the latent field parameters compressed into one x
tmb_inputs_simple_i <- tmb_inputs_simple
x_lengths <- lengths(tmb_inputs_simple_i$par_init[unique(xnames)])
tmb_inputs_simple_i$par_init$x_minus_i <- rep(0, sum(x_lengths) - 1)
tmb_inputs_simple_i$par_init$x_i <- 0
tmb_inputs_simple_i$par_init[unique(xnames)] <- NULL

#' Flag that I still have
#' $u_rho_xa
#' numeric(0)
#' in here. Not totally sure what's going on with it. I think it may not matter

#' Pass in the following as data
#' * `x_lengths`: length of each subvector within x
#' * `x_starts`: start index (remember zero indexing) of each subvector within x
#' * `i`:
x_starts <- cumsum(x_lengths) - x_lengths
tmb_inputs_simple_i$data$x_lengths <- as.numeric(x_lengths)
tmb_inputs_simple_i$data$x_starts <- as.numeric(x_starts)
tmb_inputs_simple_i$data$i <- 1

#' Change the DLL to the x_index version
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

#' Break glass in case of R crashes
# testrun <- parallel::mcparallel({})
# obj <- parallel::mccollect(testrun, wait = TRUE, timeout = 0, intermediate = FALSE)

obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = DLL,
  silent = !inner_verbose,
  random = "x_minus_i",
  map = map,
  random.start = expression(last.par[random])
)

if (!is.null(progress)) {
  obj$fn <- naomi:::report_progress(obj$fn, progress)
}

#' The set of input values that we'd like to calculate the log-probability at
#' For k = 1 there is only one row in modesandhessians so we can just pass that
(nodes <- spline_nodes(modeandhessian = modesandhessians, tmb_inputs_simple_i$data$i, k = 7))

#' Can be quite simple because there are no weights and it's just EB
laplace_marginal <- function(x, i) {
  random <- obj$env$random
  theta <- as.numeric(modesandhessians[thetanames])
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

plot_marginal_spline(nodes, lps, lower = min(nodes) - 1, upper = max(nodes) + 1)

df <- data.frame(
  index = i,
  par = xnames[i],
  x = spline_nodes(modeandhessian = modesandhessians, i = i, k = 5)
) %>%
  dplyr::mutate(lp = purrr::map_dbl(x, laplace_marginal, i = i))

#' Loop this over index i
laplace_marginal_df <- function(i) {

  data$i <- i

  obj <- TMB::MakeADFun(
    data = data,
    parameters = par,
    DLL = DLL,
    silent = !inner_verbose,
    random = "x_minus_i",
    map = map,
    random.start = expression(last.par[random])
  )

  random <- obj$env$random
  theta <- as.numeric(modesandhessians[thetanames])
  mode <- modesandhessians[["mode"]][[1]][-i]
  obj$env$last.par[random] <- mode

  f <- function(x) {
    lp <- as.numeric(- obj$fn(c(x, theta))) + log(modesandhessians$weights)
    lp - quad$normalized_posterior$lognormconst
  }

  data.frame(
    index = i,
    par = xnames[i],
    x = spline_nodes(modeandhessian = modesandhessians, i = i, k = 5)
  ) %>%
    dplyr::mutate(lp = purrr::map_dbl(x, f))
}

#' Version with DATA_UPDATE used
# laplace_marginal_df <- function(i) {
#
#   data$i <- i
#   obj$env$data$i <- i
#
#   random <- obj$env$random
#   theta <- as.numeric(modesandhessians[thetanames])
#   mode <- modesandhessians[["mode"]][[1]][-i]
#   obj$env$last.par[random] <- mode
#
#   f <- function(x) {
#     lp <- as.numeric(- obj$fn(c(x, theta))) + log(modesandhessians$weights)
#     lp - quad$normalized_posterior$lognormconst
#   }
#
#   data.frame(
#     index = i,
#     par = xnames[i],
#     x = spline_nodes(modeandhessian = modesandhessians, i = i, k = 5)
#   ) %>%
#     dplyr::mutate(lp = purrr::map_dbl(x, f))
# }

start <- Sys.time()
marginals <- purrr::map(.x = 1:10, .f = laplace_marginal_df)
end <- Sys.time()
end - start

marginals <- bind_rows(marginals)
marginals

#' How to sample from these marginals?
compute_pdf_and_cdf <- function(nodes, lps) {
  rn <- range(nodes)
  rnl <- diff(rn)
  min <- min(rn) - rnl / 2
  max <- max(rn) + rnl / 2

  ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
  interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }

  finegrid <- seq(min, max, length.out = 1000)

  df <- data.frame(
    x = finegrid,
    pdf = exp(interpolant(finegrid)),
    cdf = cumsum(exp(interpolant(finegrid))) * c(0, diff(finegrid))
  )

  return(df)
}

df <- compute_pdf_and_cdf(nodes, lps)

sample_marginal <- function(df, M) {
  M <- 1000
  q <- stats::runif(M)
  samples <- numeric(M)
  for(i in 1:M) samples[i] <- df$x[max(which(df$cdf < q[i]))]
  return(samples)
}

samples <- sample_marginal(df, M = 1000)
samples
