#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("naomi-simple")
# setwd("src/naomi-simple")

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
#' with devtools::install_github("mrc-ide/naomi", ref = "e9de40f")

compile("naomi.cpp")
dyn.load(dynlib("naomi"))

#' Want to replace the above with
# compile("naomi_simple.cpp")
# dyn.load(dynlib("naomi_simple"))

#' The replace all of the following fitting, uncertainty generation, and model output
#' with equivalent versions using naomi_simple rather than naomi. The end state should
#' be a plot of HIV prevalence 15-49 which is equivalent to the one created using naomi,
#' assuming that removing projection does not change the estimates at the time of survey
#' (and that indicators below corresponds to T1 not T2 or T3 say).

#' Fit TMB model
fit <- local_fit_tmb(tmb_inputs, outer_verbose = TRUE, inner_verbose = FALSE, max_iter = 250, progress = NULL)

#' Add uncertainty
fit <- local_sample_tmb(fit, nsample = 2000)

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
