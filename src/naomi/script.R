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
eb_quad <- sample_aghq(eb_quad, M = 2000)

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

#' 1. Four chains of 100 with four cores takes ~2.5 minutes
#' 2. Four chains of 1000 with four cores takes ~30 minutes
#' 3. Four chains of 4000 with four cores takes ~3 hours
#' I have saved the results of (3.) under the name mcmc.rds for access without waiting 3 hours!

# start_mcmc <- Sys.time()
# mcmc <- fit_tmbstan(tmb_inputs, chains = 4, iter = 4000, cores = 4)
# end_mcmc <- Sys.time()
# time_mcmc <- end_mcmc - start_mcmc

# saveRDS(mcmc, "mcmc.rds")

mcmc <- readRDS("mcmc.rds")

#' MCMC diagnostic checks

bayesplot::color_scheme_set("viridisA")
ggplot2::theme_set(theme_minimal())

#' Rhat
#' Looking for values less than 1.1 here

rhats <- bayesplot::rhat(mcmc)
bayesplot::mcmc_rhat(rhats)

#' ESS ratio
#' Worry about values less than 0.1 here

ratios <- bayesplot::neff_ratio(mcmc)
bayesplot::mcmc_neff(ratios)

#' Autocorrelation
#' High values of autocorrelation in the chains

bayesplot::mcmc_acf(mcmc, pars = vars(starts_with("beta")))

#' Univariate traceplots
#' Assess these visually

bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("beta")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("logit")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("log_sigma")))

#' Prevalence model
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_rho_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_rho_xs[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("us_rho_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("us_rho_xs[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_rho_a[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_rho_as[")))

#' ART model
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_xs[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("us_alpha_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("us_alpha_xs[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_a[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_as[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_xt[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_xa[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("u_alpha_xat[")))

#' Other
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_lambda_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_asfr_x[")))

#' ANC model
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_anc_rho_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_anc_alpha_x[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_anc_rho_xt[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("ui_anc_alpha_xt[")))
bayesplot::mcmc_trace(mcmc, pars = vars(starts_with("log_or_gamma["))) #' N.B. these are from the ANC attendance model

#' Pairs plots

#' Have suspicion that the ART attendance model is unidentifiable. Let's have a
#' look at the pairs plot for neighbouring districts and the log_or_gamma parameter.
nb <- area_merged %>%
  filter(area_level == max(area_level)) %>%
  bsae::sf_to_nb()

neighbours_log_or_gamma_pairs_plot <- function(i) {
  neighbour_pars <- paste0("log_or_gamma[", c(i, nb[[i]]), "]")
  bayesplot::mcmc_pairs(mcmc, pars = neighbour_pars, diag_fun = "hist", off_diag_fun = "hex")
}

area_merged %>%
  filter(area_level == max(area_level)) %>%
  print(n = Inf)

neighbours_log_or_gamma_pairs_plot(5)

#' NUTS specific assessment

np <- bayesplot::nuts_params(mcmc)

#' Number of divergent transitions
#' In this instance it's zero, so no need to do further divergent transition investigation
np %>%
  filter(Parameter == "divergent__") %>%
  summarise(n_divergent = sum(Value))

bayesplot::mcmc_nuts_divergence(np, bayesplot::log_posterior(mcmc))

#' Energy plots

#' Ideally these two histograms would be the same (Betancourt 2017)
#' The histograms are quite different, in a way suggesting the chains may not be
#' fully exploring the tails of the target distribution
bayesplot::mcmc_nuts_energy(np)

#' Comparison

#' All of the possible parameter names as follows
names(fit$obj$env$par)

#' Univariate plots

histogram_plot <- function(par) {
  df_compare <- rbind(
    data.frame(method = "TMB", samples = as.numeric(fit$sample[[par]])),
    data.frame(method = "aghq", samples = as.numeric(eb_quad$sample[[par]])),
    data.frame(method = "tmbstan", samples = as.numeric(unlist(rstan::extract(mcmc, pars = par))))
  )

  df_compare %>%
    group_by(method) %>%
    summarise(n = n())

  ggplot(df_compare, aes(x = samples, fill = method, col = method)) +
    geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
    theme_minimal() +
    facet_grid(method~.) +
    labs(x = paste0(par), y = "Density", fill = "Method") +
    scale_color_manual(values = multi.utils::cbpalette()) +
    scale_fill_manual(values = multi.utils::cbpalette()) +
    theme(legend.position = "none")
}

histogram_plot("beta_anc_rho")

#' Kolmogorov Smirnov tests

to_ks_df <- function(par, fit, quad, mcmc) {
  tmb <- t(fit$sample[[par]])
  aghq <- t(quad$sample[[par]])
  tmbstan <- rstan::extract(mcmc, pars = par)[[par]]

  n <- ncol(tmbstan)
  ks_tmb <- numeric(n)
  ks_aghq <- numeric(n)
  for(i in 1:n) {
    ks_tmb[i] <- inf.utils::ks_test(tmb[, i], tmbstan[, i])$D
    ks_aghq[i] <- inf.utils::ks_test(aghq[, i], tmbstan[, i])$D
  }

  rbind(
    data.frame(method = "TMB", ks = ks_tmb, par = par, index = 1:length(ks_tmb)),
    data.frame(method = "aghq", ks = ks_aghq, par = par, index = 1:length(ks_aghq))
  )
}

plot_ks_df <- function(ks_df) {
  wide_ks_df <- pivot_wider(ks_df, names_from = "method", values_from = "ks") %>%
    mutate(ks_diff = TMB - aghq)

  mean_ks_diff <- mean(wide_ks_df$ks_diff)

  boxplot <- wide_ks_df %>%
    ggplot(aes(x = ks_diff)) +
    geom_boxplot(width = 0.5) +
    coord_flip() +
    labs(
      title = paste0("Mean KS difference is ", mean_ks_diff),
      subtitle = ">0 then TMB more different to tmbstan, <0 then aghq more different",
      x = "KS(TMB, tmbstan) - KS(aghq, tmbstan)"
    ) +
    theme_minimal()

  scatterplot <- ggplot(wide_ks_df, aes(x = TMB, y = aghq)) +
    geom_point(alpha = 0.5) +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = paste0("KS tests for ", ks_df$par, " of length ", max(ks_df$index)),
      subtitle = "Values along y = x have similar KS",
      x = "KS(aghq, tmbstan)", y = "KS(TMB, tmbstan)"
    ) +
    theme_minimal()

  cowplot::plot_grid(scatterplot, boxplot, ncol = 2, rel_widths = c(1.3, 1))
}

ks_helper <- function(par) to_ks_df(par, fit = fit, quad = eb_quad, mcmc = mcmc) %>% plot_ks_df()

ks_helper("beta") #' Broken
ks_helper("logit") #' Broken
ks_helper("log_sigma") #' Broken

ks_helper("u_rho_x")
ks_helper("u_rho_xs")
ks_helper("us_rho_x")
ks_helper("us_rho_xs")
ks_helper("u_rho_a")
ks_helper("u_rho_as")

ks_helper("u_alpha_x")
ks_helper("u_alpha_xs")
ks_helper("us_alpha_x")
ks_helper("us_alpha_xs")
ks_helper("u_alpha_a")
ks_helper("u_alpha_as")
ks_helper("u_alpha_xt")
ks_helper("u_alpha_xa")
ks_helper("u_alpha_xat")

ks_helper("ui_anc_rho_x")
ks_helper("ui_anc_alpha_x")
ks_helper("ui_anc_rho_xt")
ks_helper("ui_anc_alpha_xt")
ks_helper("log_or_gamma")
