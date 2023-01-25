#' Start of if(adam) loop

#' TODO list
#' * [x] Expose internals of `fit_aghq` and run through
#' * [ ] Expose internals of `sample_aghq` to get a feel again for how the latent field and hyperparameter marginals are being used and run through
#'   * Blocker: can't fit `aghq` with `k > 1`.
#'     * Possible solution: lower `area_level` to 3 rather than 4. Won't fix `k^{#dim hyper}` problem
#'     * Possible solution: create version of `naomi_simple.cpp` with fewer `{#dim hyper}`
#' * [x] Create TMB template and objective function for one named latent field parameter and index
#' * [ ] Attain Laplace marginal for one named latent field parameter with empirical Bayes hypers
#' * [ ] Attain Laplace marginal for one named latent field parameter with dense grid hypers
#' * [ ] Attain Laplace marginal for one named latent field parameter with spare grid hypers
#'   * This may be literally the same thing as with dense grid hypers
#' * [ ] Generalise code to sweep over all indices of one named latent field parameter
#' * [ ] Generalise code to sweep over all latent field parameters
#' * [ ] Create function `fit_name_here` which brings this all together
#' * [ ] Fit `k = 1` with Laplace marginals
#' * [ ] Fit `k = 2` sparse with Laplace marginals
#' * [ ] Insert results of `k = 3` sparse into `fit_name_here` to get Laplace marginals

#' Start with EB
k <- 1

#' Start expose fit_aghq(tmb_inputs, k = k)
tmb_input <- tmb_inputs
k <- k
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
basegrid <- NULL
control <- aghq::default_control_tmb()

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

if (control$onlynormconst) return(quad)

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

quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")
quad$obj <- obj

#' modesandhessians has 1 column for each hyper, and then 1 for the mode of the latents and 1 for the Hessian of the latents
stopifnot(ncol(modesandhessians) == length(thetanames) + 2)

#' Let's begin by just trying to get the Laplace marginals working for one parameter
#' Baby steps!

compile("naomi_simple_beta_rho_index.cpp")
dyn.load(dynlib("naomi_simple_beta_rho_index"))

#' I have adapted naomi_simple_beta_rho_index to have DATA_INTEGER(i) which will
#' tell the template which index of i we are going to be Laplace approximating
#' So to do this I need to add i to the data passed in
str(tmb_inputs_simple)

#' Just put it to the first index for now
tmb_inputs_simple$data$i <- 1

#' Change the initialisation values to remove beta_rho
tmb_inputs_simple$par_init$beta_rho <- NULL
tmb_inputs_simple$par_init$beta_rho_minus_i <- 0
tmb_inputs_simple$par_init$beta_rho_i <- 0

#' Expose local_make_tmb_obj
data <- tmb_inputs_simple$data
par <- tmb_inputs_simple$par_init
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
  "beta_rho_minus_i",
  "beta_alpha",
  "beta_lambda",
  "beta_anc_rho",
  "beta_anc_alpha",
  "u_rho_x",
  "us_rho_x",
  "u_rho_xs",
  "us_rho_xs",
  "u_rho_a",
  "u_rho_as",
  "u_rho_xa",
  "u_alpha_x",
  "us_alpha_x",
  "u_alpha_xs",
  "us_alpha_xs",
  "u_alpha_a",
  "u_alpha_as",
  "u_alpha_xa",
  "ui_lambda_x",
  "ui_anc_rho_x",
  "ui_anc_alpha_x",
  "log_or_gamma"
)

obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = DLL,
  silent = !inner_verbose,
  random = integrate_out,
  map = map
)

if (!is.null(progress)) {
  obj$fn <- naomi:::report_progress(obj$fn, progress)
}
