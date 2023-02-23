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
