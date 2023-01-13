#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_fit", parameters = list(aghq = TRUE))
# setwd("src/naomi-simple_fit")

if(tmb + aghq + tmbstan != 1) {
  stop("Choose exactly one of TMB, ahgq, or tmbstan to run. Don't be greedy. I'll know.")
}

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

#' naomi.cpp was obtained from https://github.com/mrc-ide/naomi on the 7/12/22.
#' This corresponds to package number 2.8.5. You can check package version for
#' Naomi with packageVersion("naomi"), and if required, install version 2.8.5.
#' with devtools::install_github("mrc-ide/naomi", ref = "e9de40f"). We modify
#' this model to obtain a simplified version where the time points T2 and T3,
#' as well as extraneous outputs, are removed.

compile("naomi_simple.cpp")
dyn.load(dynlib("naomi_simple"))

tmb_inputs <- prepare_tmb_inputs(naomi_data)
tmb_inputs_simple <- local_exclude_inputs(tmb_inputs)

if(tmb) {
  start <- Sys.time()

  #' Fit TMB model
  fit <- local_fit_tmb(tmb_inputs_simple, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL, DLL = "naomi_simple")

  #' Add uncertainty
  fit <- local_sample_tmb(fit, nsample = nsample)

  #' Calculate model outputs
  outputs <- local_output_package_naomi_simple(fit, naomi_data)

  end <- Sys.time()

  out <- list(fit = fit, outputs = outputs, time = end - start)
  saveRDS(out, "out.rds")
}

if(aghq) {
  start <- Sys.time()

  #' The number of hyperparameters is 24, as compared with 31 for the full model
  n_hyper <- 24

  #' k = 1 empirical Bayes approach, takes ~1 minute

  #' k = 2 and ndConstruction = "sparse" it's 49 points, and takes ~45 minutes
  (mvQuad::size(mvQuad::createNIGrid(n_hyper, "GHe", 2, "sparse"))$gridpoints)

  #' k = 3 and ndConstruction = "sparse" it's 1225 points
  (mvQuad::size(mvQuad::createNIGrid(n_hyper, "GHe", 3, "sparse"))$gridpoints)

  #' Fit AGHQ model
  quad <- fit_aghq(tmb_inputs, k = k)

  #' Add uncertainty
  quad <- sample_aghq(quad, M = nsample)

  #' Error sampling from sparse quadratures at the moment due to negative weights
  # control <- aghq::default_control_tmb()
  # control$method_summaries <- "correct"
  # control$ndConstruction <- "sparse"
  # sparse_quad <- fit_aghq(tmb_inputs, k = 2, basegrid = sparse_grid_2, control = control)
  # sparse_quad <- sample_aghq(sparse_quad, M = nsample)

  end <- Sys.time()

  out <- list(quad = quad, time = end - start)
  saveRDS(out, "out.rds")
}

if(tmbstan) {
  start <- Sys.time()

  #' 1. Four chains of 100 with four cores takes ~1.5 minutes
  #' 2. Four chains of 1000 with four cores takes ~20 minutes
  #' 3. Four chains of 4000 with four cores takes ~1.5 hours
  #' 4. Four chains of 8000 with four cores takes ~3 hours

  mcmc <- fit_tmbstan(tmb_inputs, chains = 4, iter = niter, thin = nthin, cores = 4)

  end <- Sys.time()

  out <- list(mcmc = mcmc, time = end - start)
  saveRDS(out, "out.rds")
}
