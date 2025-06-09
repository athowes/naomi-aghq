#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple_fit", parameters = list(aghq = TRUE, area_level = 4, k = 3, s = 7))
# setwd("src/naomi-simple_fit")

if(tmb + aghq + adam + tmbstan != 1) {
  stop("Choose exactly one of TMB, ahgq, adam, or tmbstan to run. Don't be greedy. I'll know.")
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
level <- area_level #' Default is 4
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
#' as well as extraneous outputs, are removed. No new inference is conducted in
#' the removed model components, the difficulty of the problem remains similar.

#' This is the simplfiied version of naomi.cpp
compile("naomi_simple.cpp")
dyn.load(dynlib("naomi_simple"))

#' This is an altered version of the above where all of the latent field elements
#' are concatenated into a single vector x. It's particularly useful to do this
#' for the approach to Laplace marginals on x_i.
compile("naomi_simple_x_index.cpp")
dyn.load(dynlib("naomi_simple_x_index"))

tmb_inputs <- prepare_tmb_inputs(naomi_data)
tmb_inputs_simple <- local_exclude_inputs(tmb_inputs)

#' The number of hyperparameters is 24 (as compared with 31 for the full model)
n_hyper <- 24

#' Create version of the objective function with no Laplace approximations
#' This will be used in later reports, such as to do PSIS
objfull <- local_make_tmb_obj(tmb_inputs$data, tmb_inputs$par_init, calc_outputs = FALSE, inner_verbose = FALSE, DLL = "naomi_simple", laplace = FALSE)
saveRDS(objfull, file = "objfull.rds")

if(tmb) {
  start <- Sys.time()

  #' Fit TMB model
  fit <- local_fit_tmb(tmb_inputs_simple, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL, DLL = "naomi_simple")

  if(sample) {
    #' Add uncertainty
    fit <- local_sample_tmb(fit, random_only = random_only, M = nsample)

    #' Calculate model outputs
    outputs <- local_output_package_naomi_simple(fit, naomi_data)
  }

  end <- Sys.time()

  if(sample) {
    out <- list(fit = fit, inputs = tmb_inputs_simple, outputs = outputs, naomi_data = naomi_data, time = end - start)
  } else {
    out <- list(fit = fit, inputs = tmb_inputs_simple, naomi_data = naomi_data, time = end - start)
  }

  saveRDS(out, "out.rds")
}

if(aghq) {
  start <- Sys.time()

  inner_verbose <- FALSE
  progress <- NULL
  map <- NULL
  DLL <- "naomi_simple"

  stopifnot(inherits(tmb_inputs_simple, "naomi_simple_tmb_input"))

  obj <- local_make_tmb_obj(tmb_inputs_simple$data, tmb_inputs_simple$par_init, calc_outputs = 0L, inner_verbose, progress, map, DLL = DLL)

  #' Start by doing the optimisation for the hyperparameters
  control <- aghq::default_control_tmb()
  optresults <- optimize_theta(obj, startingvalue = obj$par, control = control)

  # In future, start optimisations from the hyperparamter mode
  startingvalue <- optresults$mode
  d <- length(optresults$mode)

  #' Now the AGHQ grid can be created, potentially using the optimisation
  stopifnot(grid_type %in% c("product", "sparse", "pca", "scaled_pca"))

  if(grid_type == "product") {
    basegrid <- NULL
    adapt <- TRUE
  }

  if(grid_type == "sparse") {
    basegrid <- mvQuad::createNIGrid(d, "GHe", k, "sparse")
    control$method_summaries <- "correct"
    control$ndConstruction <- "sparse"
    adapt <- TRUE
  }

  if(grid_type == "pca") {
    whichfirst <- 1
    idxorder <- c(whichfirst, (1:d)[-whichfirst])
    H <- optresults$hessian[idxorder, idxorder]
    basegrid <- create_pca_grid(optresults$mode[idxorder], Matrix::forceSymmetric(solve(H)), s, k)
    stopifnot(nrow(basegrid$nodes) == k^s)
    adapt <- FALSE
  }

  if(grid_type == "scaled_pca") {
    whichfirst <- 1
    idxorder <- c(whichfirst, (1:d)[-whichfirst])
    H <- optresults$hessian[idxorder, idxorder]
    basegrid <- create_scaled_pca_grid(optresults$mode[idxorder], Matrix::forceSymmetric(solve(H)), s, k)
    stopifnot(nrow(basegrid$nodes) == k^s)
    adapt <- FALSE
  }

  #' Here is where the AGHQ is done
  quad <- local_marginal_laplace_tmb(obj, k = k, optresults = optresults, startingvalue = startingvalue, basegrid = basegrid, adapt = adapt, control = control)
  objout <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 1L, inner_verbose, progress, map, DLL = DLL)
  quad$obj <- objout

  #' Add uncertainty
  #' Note: can't sample when the quadrature is sparse because of negative weights
  if(sample && grid_type != "sparse") {
    quad <- sample_aghq(quad, M = nsample)
  }

  end <- Sys.time()

  #' Note that local_output_package_naomi_simple needs to be adapted to work with aghq fits
  out <- list(quad = quad, inputs = tmb_inputs_simple, naomi_data = naomi_data, time = end - start)

  saveRDS(out, "out.rds")
}

if(adam) {
  start <- Sys.time()

  #' k = 1 empirical Bayes with Laplace marginals
  base_grid <- mvQuad::createNIGrid(dim = n_hyper, type = "GHe", level = 1, ndConstruction = "product")
  #' The k argument here shouldn't be doing anything: this should (in future) be fixed in aghq::aghq
  adam <- fit_adam(tmb_inputs_simple, k = 1, basegrid = base_grid)

  if(sample) {
    #' Add uncertainty
    adam <- sample_adam(adam, M = nsample)
  }

  #' Note that local_output_package_naomi_simple would need to be adapted to work with adam fits

  end <- Sys.time()

  out <- list(adam = adam, inputs = tmb_inputs_simple, naomi_data = naomi_data, time = end - start)

  saveRDS(out, "out.rds")
}

if(tmbstan) {
  start <- Sys.time()

  #' 1. Four chains of 100 with four cores takes ~1.5 minutes
  #' 2. Four chains of 1000 with four cores takes ~20 minutes
  #' 3. Four chains of 4000 with four cores takes ~1.5 hours
  #' 4. Four chains of 8000 with four cores takes ~3 hours

  #' Fit Stan model
  mcmc <- fit_tmbstan(tmb_inputs_simple, chains = nchains, iter = niter, thin = nthin, cores = ncores, DLL = "naomi_simple", laplace = hmc_laplace)

  #' Add uncertainty (really this is about sampling from the indicators, a.k.a. generated quantities)
  #' No M is provided here, number of samples equal to length of Markov chain are created
  #' If required, number of samples can easily be reduced afterwards
  mcmc <- sample_tmbstan(mcmc, verbose = TRUE)

  #' Note that local_output_package_naomi_simple would need to be adapted to work with tmbstan fits

  end <- Sys.time()

  out <- list(mcmc = mcmc, inputs = tmb_inputs_simple, naomi_data = naomi_data, time = end - start)
  saveRDS(out, "out.rds")
}
