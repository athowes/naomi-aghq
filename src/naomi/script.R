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
spec <- extract_pjnz_naomi(pjnz)

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

#' Expose naomi::fit_tmb
tmb_input <- tmb_inputs
outer_verbose <- TRUE
inner_verbose <- FALSE
max_iter <- 250
progress <- NULL

stopifnot(inherits(tmb_input, "naomi_tmb_input"))

compile("naomi.cpp")
dyn.load(dynlib("naomi"))

#' Expose naomi:::make_tmb_obj
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

#' Add uncertainty
fit <- sample_tmb(fit)

#' Calculate model outputs
outputs <- output_package(fit, naomi_data)

indicators <- add_output_labels(outputs) %>%
  left_join(outputs$meta_area %>% select(area_level, area_id, center_x, center_y)) %>%
  sf::st_as_sf()

pdf("15-19-prev.pdf", h = 5, w = 6.25)

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
