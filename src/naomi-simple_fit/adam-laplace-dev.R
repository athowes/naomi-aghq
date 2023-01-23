#' Start of if(adam) loop

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
