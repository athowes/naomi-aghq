#' Uncomment and run the two line below to resume development of this script
# orderly::orderly_develop_start("docs_model-structure")
# setwd("src/docs_model-structure")

tmb <- readRDS("depends/tmb.rds")
fit <- tmb$fit

#' How large is the incidence component of the model?
#' Specifically: how many latent field parameters? How many hyperparameters?

unique(names(fit$obj$env$par))

incidence_related <- c(
  "beta_lambda",
  "OmegaT_raw",
  "log_betaT",
  "log_sigma_lambda_x",
  "ui_lambda_x"
)

integrate_out <- c(
  "beta_rho", "beta_alpha", "beta_alpha_t2", "beta_lambda", "beta_asfr",
  "beta_anc_rho", "beta_anc_alpha", "beta_anc_rho_t2", "beta_anc_alpha_t2",
  "u_rho_x", "us_rho_x", "u_rho_xs", "us_rho_xs", "u_rho_a", "u_rho_as",
  "u_rho_xa", "u_alpha_x", "us_alpha_x", "u_alpha_xs", "us_alpha_xs",
  "u_alpha_a", "u_alpha_as", "u_alpha_xt", "u_alpha_xa", "u_alpha_xat",
  "u_alpha_xst", "ui_lambda_x", "logit_nu_raw", "ui_asfr_x", "ui_anc_rho_x",
  "ui_anc_alpha_x", "ui_anc_rho_xt", "ui_anc_alpha_xt", "log_or_gamma",
  "log_or_gamma_t1t2"
)

exclude_random_pars <- c(
  "beta_alpha_t2", "beta_asfr", "beta_anc_rho_t2", "beta_anc_alpha_t2",
  "u_alpha_xt", "u_alpha_xat", "u_alpha_xst", "logit_nu_raw", "ui_asfr_x",
  "ui_anc_rho_xt", "ui_anc_alpha_xt", "log_or_gamma_t1t2"
)

integrate_out <- setdiff(integrate_out, exclude_random_pars)
latent_field <- integrate_out

hyper <- setdiff(unique(names(fit$obj$env$par)), latent_field)

#' 3 out of 24 hypers are in the incidence component
sum(incidence_related %in% hyper)
length(hyper)

#' 34 out of 467 latent field are in the incidence component
sum(names(fit$obj$env$par) %in% incidence_related[incidence_related %in% latent_field])
length(names(fit$obj$env$par)) - length(hyper)

#' Analysis of inputs to Naomi outputs

#' Inputs for rho_t1_out
# vector<Type> rho_t1_out(plhiv_t1_out / population_t1_out);
# vector<Type> plhiv_t1_out(A_out * plhiv_t1)
# vector<Type> plhiv_t1(population_t1 * rho_t1);

#' Inputs for alpha_t1_out
# vector<Type> alpha_t1_out(artnum_t1_out / plhiv_t1_out);
# vector<Type> artnum_t1_out(A_out * artnum_t1);
# vector<Type> artnum_t1(population_t1 * prop_art_t1);

#' Inputs for lambda_t1_out
# vector<Type> lambda_t1_out(infections_t1_out / (population_t1_out - plhiv_t1_out));
# vector<Type> infections_t1_out(A_out * infections_t1);
# vector<Type> infections_t1(lambda_t1 * (population_t1 - plhiv_t1));

# Add to Table S2?
