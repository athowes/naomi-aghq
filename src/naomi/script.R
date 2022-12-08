#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi")
# setwd("src/naomi")

#' Obtain data
area_merged <- read_sf(system.file("extdata/demo_areas.geojson", package = "naomi"))
pop_agesex <- read_csv(system.file("extdata/demo_population_agesex.csv", package = "naomi"))
survey_hiv_indicators <- read_csv(system.file("extdata/demo_survey_hiv_indicators.csv", package = "naomi"))
art_number <- read_csv(system.file("extdata/demo_art_number.csv", package = "naomi"))
anc_testing <- read_csv(system.file("extdata/demo_anc_testing.csv", package = "naomi"))
pjnz <- system.file("extdata/demo_mwi2019.PJNZ", package = "naomi")
spec <- naomi::extract_pjnz_naomi(pjnz)

#' Settings
scope <- "MWI"
level <- 4
calendar_quarter_t1 <- "CY2016Q1"
calendar_quarter_t2 <- "CY2018Q3"
calendar_quarter_t3 <- "CY2019Q4"
prev_survey_ids  <- c("DEMO2016PHIA", "DEMO2015DHS")
artcov_survey_ids  <- "DEMO2016PHIA"
vls_survey_ids <- NULL
recent_survey_ids <- "DEMO2016PHIA"
artnum_calendar_quarter_t1 <- "CY2016Q1"
artnum_calendar_quarter_t2 <- "CY2018Q3"
anc_clients_year2 <- 2018
anc_clients_year2_num_months <- 9
anc_prevalence_year1 <- 2016
anc_prevalence_year2 <- 2018
anc_art_coverage_year1 <- 2016
anc_art_coverage_year2 <- 2018

#' Prepare model inputs
naomi_mf <- naomi_model_frame(
  area_merged,
  pop_agesex,
  spec,
  scope = scope,
  level = level,
  calendar_quarter_t1,
  calendar_quarter_t2,
  calendar_quarter_t3
)

naomi_data <- select_naomi_data(
  naomi_mf,
  survey_hiv_indicators,
  anc_testing,
  art_number,
  prev_survey_ids,
  artcov_survey_ids,
  recent_survey_ids,
  vls_survey_ids,
  artnum_calendar_quarter_t1,
  artnum_calendar_quarter_t2,
  anc_prevalence_year1,
  anc_prevalence_year2,
  anc_art_coverage_year1,
  anc_art_coverage_year2
)

#' Fit model with TMB
tmb_inputs <- prepare_tmb_inputs(naomi_data)

#' Note that naomi.cpp was obtained from https://github.com/mrc-ide/naomi on the 7/12/22
#' This corresponds to package number 2.8.5.

#' Check package version for Naomi -- you probably want it to match 2.8.5.
packageVersion("naomi")

#' Version 2.8.5. of Naomi can be installed with (TODO: how to select version number)
# devtools::install_github("mrc-ide/naomi")

compile("naomi.cpp")
dyn.load(dynlib("naomi"))

#' Begin expose naomi::fit_tmb
#' https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L557
#' Replacing the following call
#' fit <- fit_tmb(tmb_inputs, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL)
tmb_input <- tmb_inputs
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL

stopifnot(inherits(tmb_input, "naomi_tmb_input"))

#' Begin expose naomi:::make_tmb_obj
#' https://github.com/mrc-ide/naomi/blob/e9de40f12cf2e652f78966bb351fa5718ecd7867/R/tmb-model.R#L496
#' Replacing the following call
#' obj <- naomi:::make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress)
data <- tmb_input$data
par <- tmb_input$par_init
calc_outputs <- 0L

data$calc_outputs <- as.integer(calc_outputs)

obj <- TMB::MakeADFun(
  data = data,
  parameters = par,
  DLL = "naomi",
  silent = !inner_verbose,
  random = c("beta_rho", "beta_alpha",
    "beta_alpha_t2", "beta_lambda", "beta_asfr", "beta_anc_rho",
    "beta_anc_alpha", "beta_anc_rho_t2", "beta_anc_alpha_t2",
    "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs", "u_rho_a",
    "u_rho_as", "u_rho_xa", "u_alpha_x", "us_alpha_x",
    "u_alpha_xs", "us_alpha_xs", "u_alpha_a", "u_alpha_as",
    "u_alpha_xt", "u_alpha_xa", "u_alpha_xat", "u_alpha_xst",
    "ui_lambda_x", "logit_nu_raw", "ui_asfr_x", "ui_anc_rho_x",
    "ui_anc_alpha_x", "ui_anc_rho_xt", "ui_anc_alpha_xt",
    "log_or_gamma", "log_or_gamma_t1t2")
)

if (!is.null(progress)) {
  obj$fn <- report_progress(obj$fn, progress)
}

#' End expose naomi:::make_tmb_obj

trace <- ifelse(outer_verbose, 1, 0)

f <- withCallingHandlers(
  stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace, iter.max = max_iter)),
  warning = function(w) {
    if (grepl("NA/NaN function evaluation", w$message)) invokeRestart("muffleWarning")
  }
)

if (f$convergence != 0) warning(paste("convergence error:", f$message))

if (outer_verbose) message(paste("converged:", f$message))

f$par.fixed <- f$par
f$par.full <- obj$env$last.par
objout <- naomi:::make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 1L, inner_verbose)
f$mode <- objout$report(f$par.full)
val <- c(f, obj = list(objout))
class(val) <- "naomi_fit"

fit <- val

#' End expose naomi::fit_tmb

#' Add uncertainty
fit <- sample_tmb(fit)

#' Calculate model outputs
outputs <- output_package(fit, naomi_data)

indicators <- add_output_labels(outputs) %>%
  left_join(outputs$meta_area %>% select(area_level, area_id, center_x, center_y)) %>%
  sf::st_as_sf()

pdf("tmb-15-19-prev.pdf", h = 5, w = 6.25)

indicators %>%
  filter(age_group == "Y015_049",
         indicator == "prevalence",
         area_level == 4) %>%
  ggplot(aes(fill = mode)) +
  geom_sf() +
  viridis::scale_fill_viridis(labels = scales::percent_format()) +
  th_map() +
  facet_wrap(~sex)

dev.off()

#' Develop new function naomi::fit_aghq
#' fit <- fit_aghq(tmb_inputs, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL)
tmb_input <- tmb_inputs
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL

stopifnot(inherits(tmb_input, "naomi_tmb_input"))

obj <- naomi:::make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress)

#' The number of hyperparameters is 31
length(obj$par)

#' This is quite a lot if we plan on using a dense grid
#' For example with k = 3 points per dimension then we'll get 3^31 ~= 6E14 points
#' Or with k = 2 points per dimension 2^31 ~= 2E9
#' Let's set k = 1 (empirical Bayes) for now, and get all of the infrastructure working first

#' Using the branch https://github.com/athowes/aghq/tree/issue6 for now on this to fix issue with k = 1
#' Will PR this to master soon!
quad <- aghq::marginal_laplace_tmb(obj, k = 1, startingvalue = obj$par)
