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

tmb_inputs <- prepare_tmb_inputs(naomi_data)

#' naomi.cpp was obtained from https://github.com/mrc-ide/naomi on the 7/12/22.
#' This corresponds to package number 2.8.5. You can check package version for
#' Naomi with packageVersion("naomi"), and if required, install version 2.8.5.
#' with devtools::install_github("mrc-ide/naomi", ref = "e9de40f")

compile("naomi.cpp")
dyn.load(dynlib("naomi"))

#' Fit TMB model
fit <- local_fit_tmb(tmb_inputs, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL)

#' Add uncertainty
fit <- local_sample_tmb(fit, nsample = 100)

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

#' aghq

#' The number of hyperparameters is 31, which is a lot if we plan on using a
#' dense grid. With k = 3 points per dimension then we'll get 3^31 ~= 6E14 points
#' Or with k = 2 points per dimension 2^31 ~= 2E9. Let's set k = 1 (empirical Bayes)
#' for now, and get all of the infrastructure working first
start_eb_quad <- Sys.time()
eb_quad <- fit_aghq(tmb_inputs, k = 1)
end_eb_quad <- Sys.time()
time_eb_quad <- end_eb_quad - start_eb_quad

#' Eventually want to make aghq output match that of TMB for the Naomi model.
#' Should be able to feed into the functions naomi::output_package, naomi:::extract_indicators
#' and naomi::extract_art_attendance. The output of naomi::fit_tmb is pretty complex,
#' but sample_tmb just adds $sample, a list where each component is a matrix of samples
#' from the named parameter.
str(fit)
str(fit$sample)

#' Add uncertainty
eb_quad <- sample_aghq(eb_quad, M = 100)

#' With k = 2 and ndConstruction = "sparse" it's 63 points: should be feasible
n_hyper <- length(fit$obj$env$par) - length(fit$obj$env$random)
sparse_grid <- mvQuad::createNIGrid(n_hyper, "GHe", 2, "sparse")
mvQuad::size(sparse_grid)$gridpoints

#' Set aghq only to compute the normalising constant. For sparse grids some downsteam
#' elements, like the hyperparameter marginals, are not working yet because a product
#' structure is assumed. Unsure if this is a real limitation, or just something which
#' could be implemented but hasn't yet.
control <- aghq::default_control_tmb()
control$onlynormconst <- TRUE

start_sparse_quad <- Sys.time()
sparse_quad <- fit_aghq(tmb_inputs, k = 2, basegrid = sparse_grid, control = control)
end_sparse_quad <- Sys.time()
time_sparse_quad <- end_sparse_quad - start_sparse_quad

#' tmbstan

mcmc <- fit_tmbstan(tmb_inputs, chains = 4, iter = 100, warmup = 75)

#' Comparison
names(fit$obj$env$par)
test_par <- "beta_anc_rho" #' Picked pretty randomly
df_compare <- rbind(
  data.frame(method = "TMB", samples = as.numeric(fit$sample[[test_par]])),
  data.frame(method = "aghq", samples = as.numeric(eb_quad$sample[[test_par]])),
  data.frame(method = "tmbstan", samples = as.numeric(unlist(rstan::extract(mcmc, pars = test_par))))
)

df_compare %>%
  group_by(method) %>%
  summarise(n = n())

ggplot(df_compare, aes(x = samples, fill = method)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = paste0(test_par), y = "Density", fill = "Method")
