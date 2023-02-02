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
