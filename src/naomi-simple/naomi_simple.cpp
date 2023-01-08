// #define TMB_LIB_INIT R_init_naomi_simple
#include <TMB.hpp>


/** Log posterior density of BYM2 with INLA conditional parameterisation
 *
 * Calculate the joint LPDF of parameter vector (x, u) where
 * x = sigma * (sqrt(phi) * u + sqrt(1-phi) * v) with u a ICAR structured
 * component u ~ N(0, Q^-1) and v is an IID effect v ~ N(0, 1). Calculation
 * proceeds by conditioning P(<x, u>) = P(x | u) * P(u). See Reibler et al.
 * Section 3.4.
 *
 * @param x vector of random effects.
 * @param u vector of spatial component of random effect.
 * @param sigma marginal standard deviation (>0).
 * @param phi proportion of marginal variance explained by spatial structured
 *            component u (phi \in [0, 1]).
 * @param Q scaled structure matrix for spatial component.
 *
 * @return Log probability density of x and u.
 *
 * @note
 * The $\sqrt(2\pi)^{-2*n}$ and $|Q|^{1/2}$ terms are dropped.
 * Returns the _positive_ log PDF (different from builtin TMB
 * functions. Thus should typically be implemented as `nll -= bym2_conditional_lpdf(...)`.
 */
template<class Type>
Type bym2_conditional_lpdf(const vector<Type> x,
                           const vector<Type> u,
                           const Type sigma,
                           const Type phi,
                           const Eigen::SparseMatrix<Type> Q) {

  Type val(0.0);

  // constant terms omitted: -0.5 * (n + rank(Q)) * log(2*pi) + 0.5 * log|Q|
  val += -0.5 * x.size() * (2 * log(sigma) + log(1 - phi));  // normalising constant
  val += -0.5 / (sigma * sigma * (1 - phi)) * (x * x).sum();
  val += sqrt(phi) / (sigma * (1 - phi)) * (x * u).sum();
  val += -0.5 * (u * (Q * u)).sum();
  val += -0.5 * phi / (1 - phi) * (u * u).sum();

  return(val);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // indexing:
  //
  // * rho: HIV prevalence model
  // * alpha: ART coverage model
  // * lambda: HIV incidence model
  //
  // * _x: area
  // * _a: age
  // * _s: sex
  // * _t: time

  using namespace density;

  // ** Data **

  // Population
  DATA_VECTOR(population_t1);

  // Design matrices
  DATA_MATRIX(X_rho);
  DATA_MATRIX(X_alpha);
  DATA_MATRIX(X_lambda);

  DATA_MATRIX(X_ancrho);
  DATA_MATRIX(X_ancalpha);

  DATA_SPARSE_MATRIX(Z_rho_x);
  DATA_SPARSE_MATRIX(Z_rho_xs);
  DATA_SPARSE_MATRIX(Z_rho_a);
  DATA_SPARSE_MATRIX(Z_rho_as);
  DATA_SPARSE_MATRIX(Z_rho_xa);

  DATA_SPARSE_MATRIX(Z_alpha_x);
  DATA_SPARSE_MATRIX(Z_alpha_xs);
  DATA_SPARSE_MATRIX(Z_alpha_a);
  DATA_SPARSE_MATRIX(Z_alpha_as);
  DATA_SPARSE_MATRIX(Z_alpha_xa);

  DATA_SPARSE_MATRIX(Z_x);
  DATA_SPARSE_MATRIX(Z_lambda_x);

  DATA_VECTOR(logit_rho_offset);
  DATA_VECTOR(logit_alpha_offset);

  DATA_VECTOR(log_asfr_t1_offset);

  DATA_VECTOR(logit_anc_rho_t1_offset);

  DATA_VECTOR(logit_anc_alpha_t1_offset);

  DATA_SPARSE_MATRIX(Z_ancrho_x);
  DATA_SPARSE_MATRIX(Z_ancalpha_x);


  // Precision matrix for ICAR area model
  DATA_SPARSE_MATRIX(Q_x);
  DATA_SCALAR(Q_x_rankdef);

  DATA_VECTOR(n_prev);
  DATA_VECTOR(x_prev);
  DATA_SPARSE_MATRIX(A_prev);

  DATA_VECTOR(n_artcov);
  DATA_VECTOR(x_artcov);
  DATA_SPARSE_MATRIX(A_artcov);

  DATA_VECTOR(n_recent);
  DATA_VECTOR(x_recent);
  DATA_SPARSE_MATRIX(A_recent);

  DATA_VECTOR(n_anc_prev_t1);
  DATA_VECTOR(x_anc_prev_t1);
  DATA_SPARSE_MATRIX(A_anc_prev_t1);

  DATA_VECTOR(n_anc_artcov_t1);
  DATA_VECTOR(x_anc_artcov_t1);
  DATA_SPARSE_MATRIX(A_anc_artcov_t1);

  DATA_SPARSE_MATRIX(A_artattend_t1);
  DATA_VECTOR(x_artnum_t1);

  DATA_SPARSE_MATRIX(A_artattend_mf);
  DATA_SPARSE_MATRIX(A_art_reside_attend);

  DATA_IVECTOR(n_nb);
  DATA_IVECTOR(adj_i);
  DATA_IVECTOR(adj_j);

  DATA_SPARSE_MATRIX(Xgamma);
  DATA_VECTOR(log_gamma_offset);

  DATA_SPARSE_MATRIX(Xart_idx);
  DATA_SPARSE_MATRIX(Xart_gamma);


  // Incidence model
  DATA_SCALAR(omega);
  DATA_SCALAR(OmegaT0);
  DATA_SCALAR(sigma_OmegaT);
  DATA_SCALAR(betaT0);
  DATA_SCALAR(sigma_betaT);
  DATA_SCALAR(ritaT);

  DATA_SPARSE_MATRIX(X_15to49);
  DATA_VECTOR(log_lambda_t1_offset);

  // Paediatric prevalence and incidence ratio model

  DATA_SPARSE_MATRIX(X_15to49f);
  DATA_SPARSE_MATRIX(X_paed_rho_ratio);
  DATA_VECTOR(paed_rho_ratio_offset);

  DATA_SPARSE_MATRIX(X_paed_lambda_ratio_t1);

  // ** Initialize nll **
  Type val(0);

  // ** Parameters **

  // fixed effects
  // diffuse N(0.0, 5.0) prior distribution

  PARAMETER_VECTOR(beta_rho);
  val -= dnorm(beta_rho, 0.0, 5.0, true).sum();

  PARAMETER_VECTOR(beta_alpha);
  val -= dnorm(beta_alpha, 0.0, 5.0, true).sum();

  PARAMETER_VECTOR(beta_lambda);
  val -= dnorm(beta_lambda, 0.0, 5.0, true).sum();

  PARAMETER_VECTOR(beta_anc_rho);
  val -= dnorm(beta_anc_rho, 0.0, 5.0, true).sum();

  PARAMETER_VECTOR(beta_anc_alpha);
  val -= dnorm(beta_anc_alpha, 0.0, 5.0, true).sum();

  // * HIV prevalence model *

  // hyper parameters

  PARAMETER(logit_phi_rho_x);
  Type phi_rho_x(invlogit(logit_phi_rho_x));
  val -= log(phi_rho_x) +  log(1 - phi_rho_x);  // change of variables: logit_phi_x ->v phi_x
  val -= dbeta(phi_rho_x, Type(0.5), Type(0.5), true);

  PARAMETER(log_sigma_rho_x);
  Type sigma_rho_x(exp(log_sigma_rho_x));
  val -= dnorm(sigma_rho_x, Type(0.0), Type(2.5), true) + log_sigma_rho_x;

  PARAMETER(logit_phi_rho_xs);
  Type phi_rho_xs(invlogit(logit_phi_rho_xs));
  val -= log(phi_rho_xs) +  log(1 - phi_rho_xs);  // change of variables: logit_phi_xs -> phi_xs
  val -= dbeta(phi_rho_xs, Type(0.5), Type(0.5), true);

  PARAMETER(log_sigma_rho_xs);
  Type sigma_rho_xs(exp(log_sigma_rho_xs));
  val -= dnorm(sigma_rho_xs, Type(0.0), Type(2.5), true) + log_sigma_rho_xs;

  PARAMETER(logit_phi_rho_a);
  val -= dnorm(logit_phi_rho_a, Type(0.0), Type(2.582), true);  // INLA default
  Type phi_rho_a(2.0 * invlogit(logit_phi_rho_a) - 1.0);

  PARAMETER(log_sigma_rho_a);
  Type sigma_rho_a(exp(log_sigma_rho_a));
  val -= dnorm(sigma_rho_a, Type(0.0), Type(2.5), true) + log_sigma_rho_a;

  PARAMETER(logit_phi_rho_as);
  val -= dnorm(logit_phi_rho_as, Type(0.0), Type(2.582), true);  // INLA default
  Type phi_rho_as(2.0 * invlogit(logit_phi_rho_as) - 1.0);

  PARAMETER(log_sigma_rho_as);
  Type sigma_rho_as(exp(log_sigma_rho_as));
  val -= dnorm(sigma_rho_as, Type(0.0), Type(2.5), true) + log_sigma_rho_as;

  PARAMETER(log_sigma_rho_xa);
  Type sigma_rho_xa(exp(log_sigma_rho_xa));
  val -= dnorm(sigma_rho_xa, Type(0.0), Type(0.5), true) + log_sigma_rho_xa;

  // latent effects

  PARAMETER_VECTOR(u_rho_x);
  PARAMETER_VECTOR(us_rho_x);
  val -= dnorm(sum(us_rho_x), Type(0.0), Type(0.001) * us_rho_x.size(), true); // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(u_rho_x, us_rho_x, sigma_rho_x, phi_rho_x, Q_x);


  PARAMETER_VECTOR(u_rho_xs);
  PARAMETER_VECTOR(us_rho_xs);
  if (u_rho_xs.size()) {
    val -= dnorm(sum(us_rho_xs), Type(0.0), Type(0.001) * us_rho_xs.size(), true); // soft sum-to-zero constraint
    val -= bym2_conditional_lpdf(u_rho_xs, us_rho_xs, sigma_rho_xs, phi_rho_xs, Q_x);
  }

  PARAMETER_VECTOR(u_rho_a);
  if(u_rho_a.size() > 0)
    val += SCALE(AR1(phi_rho_a), sigma_rho_a)(u_rho_a);

  PARAMETER_VECTOR(u_rho_as);
  if(u_rho_a.size() > 0)
    val += SCALE(AR1(phi_rho_as), sigma_rho_as)(u_rho_as);

  PARAMETER_VECTOR(u_rho_xa);
  if (u_rho_xa.size() > 0) {
    val -= dnorm(sum(u_rho_xa), Type(0.0), sigma_rho_xa * Type(0.001) * u_rho_xa.size(), true); // soft sum-to-zero constraint

    val -= -(Q_x.rows() - Q_x_rankdef) * log_sigma_rho_xa -
      0.5 / (sigma_rho_xa * sigma_rho_xa) * (u_rho_xa * (Q_x * u_rho_xa)).sum();
  }


  // * ART coverage model *

  PARAMETER(logit_phi_alpha_x);
  Type phi_alpha_x(invlogit(logit_phi_alpha_x));
  val -= log(phi_alpha_x) +  log(1 - phi_alpha_x);  // change of variables: logit_phi_x -> phi_x
  val -= dbeta(phi_alpha_x, Type(0.5), Type(0.5), true);

  PARAMETER(log_sigma_alpha_x);
  Type sigma_alpha_x(exp(log_sigma_alpha_x));
  val -= dnorm(sigma_alpha_x, Type(0.0), Type(2.5), true) + log_sigma_alpha_x;

  PARAMETER(logit_phi_alpha_xs);
  Type phi_alpha_xs(invlogit(logit_phi_alpha_xs));
  val -= log(phi_alpha_xs) +  log(1 - phi_alpha_xs);  // change of variables: logit_phi_xs -> phi_xs
  val -= dbeta(phi_alpha_xs, Type(0.5), Type(0.5), true);

  PARAMETER(log_sigma_alpha_xs);
  Type sigma_alpha_xs(exp(log_sigma_alpha_xs));
  val -= dnorm(sigma_alpha_xs, Type(0.0), Type(2.5), true) + log_sigma_alpha_xs;

  PARAMETER(logit_phi_alpha_a);
  val -= dnorm(logit_phi_alpha_a, Type(0.0), Type(2.582), true);  // INLA default
  Type phi_alpha_a(2.0 * invlogit(logit_phi_alpha_a) - 1.0);

  PARAMETER(log_sigma_alpha_a);
  Type sigma_alpha_a(exp(log_sigma_alpha_a));
  val -= dnorm(sigma_alpha_a, Type(0.0), Type(2.5), true) + log_sigma_alpha_a;

  PARAMETER(logit_phi_alpha_as);
  val -= dnorm(logit_phi_alpha_as, Type(0.0), Type(2.582), true);  // INLA default
  Type phi_alpha_as(2.0 * invlogit(logit_phi_alpha_as) - 1.0);

  PARAMETER(log_sigma_alpha_as);
  Type sigma_alpha_as(exp(log_sigma_alpha_as));
  val -= dnorm(sigma_alpha_as, Type(0.0), Type(2.5), true) + log_sigma_alpha_as;

  PARAMETER(log_sigma_alpha_xa);
  Type sigma_alpha_xa(exp(log_sigma_alpha_xa));
  val -= dnorm(sigma_alpha_xa, Type(0.0), Type(2.5), true) + log_sigma_alpha_xa;

  PARAMETER_VECTOR(u_alpha_x);
  PARAMETER_VECTOR(us_alpha_x);
  val -= dnorm(sum(us_alpha_x), Type(0.0), Type(0.001) * us_alpha_x.size(), true); // soft sum-to-zero constraint
  val -= bym2_conditional_lpdf(u_alpha_x, us_alpha_x, sigma_alpha_x, phi_alpha_x, Q_x);

  PARAMETER_VECTOR(u_alpha_xs);
  PARAMETER_VECTOR(us_alpha_xs);
  if (u_alpha_xs.size()) {
    val -= dnorm(sum(us_alpha_xs), Type(0.0), Type(0.001) * us_alpha_xs.size(), true); // soft sum-to-zero constraint
    val -= bym2_conditional_lpdf(u_alpha_xs, us_alpha_xs, sigma_alpha_xs, phi_alpha_xs, Q_x);
  }

  PARAMETER_VECTOR(u_alpha_a);
  if(u_alpha_a.size() > 0)
    val += SCALE(AR1(phi_alpha_a), sigma_alpha_a)(u_alpha_a);

  PARAMETER_VECTOR(u_alpha_as);
  if(u_alpha_as.size() > 0)
    val += SCALE(AR1(phi_alpha_as), sigma_alpha_as)(u_alpha_as);

  PARAMETER_VECTOR(u_alpha_xa);
  val -= dnorm(u_alpha_xa, 0.0, sigma_alpha_xa, true).sum();


  // * HIV incidence model *

  PARAMETER(OmegaT_raw);
  val -= dnorm(OmegaT_raw, Type(0.0), Type(1.0), true);
  Type OmegaT = OmegaT0 + OmegaT_raw * sigma_OmegaT;

  PARAMETER(log_betaT);
  val -= dnorm(exp(log_betaT), Type(0.0), Type(1.0), true) + log_betaT;
  Type betaT = betaT0 + exp(log_betaT) * sigma_betaT;

  PARAMETER(log_sigma_lambda_x);
  Type sigma_lambda_x(exp(log_sigma_lambda_x));
  val -= dnorm(sigma_lambda_x, Type(0.0), Type(1.0), true) + log_sigma_lambda_x;

  PARAMETER_VECTOR(ui_lambda_x);
  val -= sum(dnorm(ui_lambda_x, 0.0, sigma_lambda_x, true));


  // * ANC testing model *

  // ANC prevalence and ART coverage random effects
  PARAMETER(log_sigma_ancrho_x);
  Type sigma_ancrho_x(exp(log_sigma_ancrho_x));
  val -= dnorm(sigma_ancrho_x, Type(0.0), Type(2.5), true) + log_sigma_ancrho_x;

  PARAMETER(log_sigma_ancalpha_x);
  Type sigma_ancalpha_x(exp(log_sigma_ancalpha_x));
  val -= dnorm(sigma_ancalpha_x, Type(0.0), Type(2.5), true) + log_sigma_ancalpha_x;

  PARAMETER_VECTOR(ui_anc_rho_x);
  val -= sum(dnorm(ui_anc_rho_x, 0.0, sigma_ancrho_x, true));

  PARAMETER_VECTOR(ui_anc_alpha_x);
  val -= sum(dnorm(ui_anc_alpha_x, 0.0, sigma_ancalpha_x, true));


  // * ART attendance model *

  PARAMETER(log_sigma_or_gamma);
  Type sigma_or_gamma(exp(log_sigma_or_gamma));
  val -= dnorm(sigma_or_gamma, Type(0.0), Type(2.5), true) + log_sigma_or_gamma;

  PARAMETER_VECTOR(log_or_gamma);
  val -= dnorm(log_or_gamma, 0.0, sigma_or_gamma, true).sum();


  // *** Process model ***

  // HIV prevalence time 1

  vector<Type> mu_rho(X_rho * beta_rho +
    logit_rho_offset +
    Z_rho_x * u_rho_x +
    Z_rho_xs * u_rho_xs +
    Z_rho_a * u_rho_a +
    Z_rho_as * u_rho_as +
    Z_rho_xa * u_rho_xa);

  // paediatric prevalence

  vector<Type> rho_15to49f_t1((X_15to49f * vector<Type>(invlogit(mu_rho) * population_t1)) / (X_15to49f * population_t1));
  vector<Type> mu_rho_paed(X_paed_rho_ratio * rho_15to49f_t1 + paed_rho_ratio_offset);
  mu_rho_paed = logit(mu_rho_paed);
  mu_rho += mu_rho_paed;

  // ART coverage time 1

  vector<Type> mu_alpha(X_alpha * beta_alpha +
    logit_alpha_offset +
    Z_alpha_x * u_alpha_x +
    Z_alpha_xs * u_alpha_xs +
    Z_alpha_a * u_alpha_a +
    Z_alpha_as * u_alpha_as +
    Z_alpha_xa * u_alpha_xa);


  vector<Type> rho_t1(invlogit(mu_rho));
  vector<Type> alpha_t1(invlogit(mu_alpha));

  vector<Type> plhiv_t1(population_t1 * rho_t1);
  vector<Type> prop_art_t1(rho_t1 * alpha_t1);
  vector<Type> artnum_t1(population_t1 * prop_art_t1);

  vector<Type> plhiv_15to49_t1(X_15to49 * plhiv_t1);
  vector<Type> rho_15to49_t1(plhiv_15to49_t1 / (X_15to49 * population_t1));
  vector<Type> alpha_15to49_t1((X_15to49 * artnum_t1) / plhiv_15to49_t1);

  vector<Type> mu_lambda_t1(X_lambda * beta_lambda + log_lambda_t1_offset +
    Z_x * vector<Type>(log(rho_15to49_t1) + log(1.0 - omega * alpha_15to49_t1)) +
    Z_lambda_x * ui_lambda_x);

  vector<Type> lambda_adult_t1(exp(mu_lambda_t1));

  // Add paediatric incidence
  vector<Type> lambda_paed_t1(X_paed_lambda_ratio_t1 * rho_15to49f_t1);
  vector<Type> lambda_t1(lambda_adult_t1 + lambda_paed_t1);

  vector<Type> infections_t1(lambda_t1 * (population_t1 - plhiv_t1));

  // likelihood for household survey data

  vector<Type> rho_obs_t1((A_prev * plhiv_t1) / (A_prev * population_t1));
  vector<Type> hhs_prev_ll = dbinom(x_prev, n_prev, rho_obs_t1, true);
  val -= sum(hhs_prev_ll);

  vector<Type> alpha_obs_t1((A_artcov * artnum_t1) / (A_artcov * plhiv_t1));
  vector<Type> hhs_artcov_ll = dbinom(x_artcov, n_artcov, alpha_obs_t1, true);
  val -= sum(hhs_artcov_ll);

  vector<Type> pR_infections_obs_t1(A_recent * infections_t1);
  vector<Type> pR_plhiv_obs_t1(A_recent * plhiv_t1);
  vector<Type> pR_population_obs_t1(A_recent * population_t1);
  vector<Type> pR_lambda_obs_t1(pR_infections_obs_t1 / (pR_population_obs_t1 - pR_plhiv_obs_t1));
  vector<Type> pR_rho_obs_t1(pR_plhiv_obs_t1 / pR_population_obs_t1);
  vector<Type> pR(1.0 - exp(-(pR_lambda_obs_t1 * (1.0 - pR_rho_obs_t1) / pR_rho_obs_t1 *
    (OmegaT - betaT * ritaT) + betaT)));
  val -= dbinom(x_recent, n_recent, pR, true).sum();


  // ANC prevalence and ART coverage model
  // Note: currently this operates on the entire population vector, producing
  //       lots of zeros for males and female age groups not exposed to fertility.
  //       It would be more computationally efficient to project this to subset
  //       of female age 15-49 age groups. But I don't know if it would be
  //       meaningfully more efficient.

  vector<Type> mu_anc_rho_t1(mu_rho +
    logit_anc_rho_t1_offset +
    X_ancrho * beta_anc_rho +
    Z_ancrho_x * ui_anc_rho_x);
  vector<Type> anc_rho_t1(invlogit(mu_anc_rho_t1));

  vector<Type> mu_anc_alpha_t1(mu_alpha +
    logit_anc_alpha_t1_offset +
    X_ancalpha * beta_anc_alpha +
    Z_ancalpha_x * ui_anc_alpha_x);
  vector<Type> anc_alpha_t1(invlogit(mu_anc_alpha_t1));

  // JE NOTE 6 Jan 2022: removed mu_asfr term -- should not use for aggregate ANC.
  vector<Type> anc_clients_t1(population_t1 * exp(log_asfr_t1_offset));
  vector<Type> anc_plhiv_t1(anc_clients_t1 * anc_rho_t1);
  vector<Type> anc_already_art_t1(anc_plhiv_t1 * anc_alpha_t1);

  // likelihood for ANC testing observations

  vector<Type> anc_rho_obs_t1(A_anc_prev_t1 * anc_plhiv_t1 / (A_anc_prev_t1 * anc_clients_t1));
  vector<Type> anc_rho_obs_t1_ll = dbinom(x_anc_prev_t1, n_anc_prev_t1, anc_rho_obs_t1, true);
  val -= sum(anc_rho_obs_t1_ll);

  vector<Type> anc_alpha_obs_t1(A_anc_artcov_t1 * anc_already_art_t1 / (A_anc_artcov_t1 * anc_plhiv_t1));
  vector<Type> anc_alpha_obs_t1_ll = dbinom(x_anc_artcov_t1, n_anc_artcov_t1, anc_alpha_obs_t1, true);
  val -= sum(anc_alpha_obs_t1_ll);

  // * ART attendance model *

  vector<Type> gamma_art_t1(exp(Xgamma * log_or_gamma + log_gamma_offset));
  int cum_nb = 0;
  for(int i = 0; i < n_nb.size(); i++){
    Type cum_exp_or_gamma_i = 0.0;
    for(int j = 0; j < n_nb[i]+1; j++)
      cum_exp_or_gamma_i += gamma_art_t1[cum_nb + i + j];
    for(int j = 0; j < n_nb[i]+1; j++)
      gamma_art_t1[cum_nb + i + j] /= cum_exp_or_gamma_i;
    cum_nb += n_nb[i];
  }

  vector<Type> prop_art_ij_t1((Xart_idx * prop_art_t1) * (Xart_gamma * gamma_art_t1));
  vector<Type> population_ij_t1(Xart_idx * population_t1);

  vector<Type> artnum_ij_t1(population_ij_t1 * prop_art_ij_t1);
  vector<Type> A_j_t1(A_artattend_t1 * artnum_ij_t1);
  vector<Type> sd_A_j_t1(A_artattend_t1 * vector<Type>(population_ij_t1 * prop_art_ij_t1 * (1 - prop_art_ij_t1)));
  sd_A_j_t1 = sd_A_j_t1.sqrt();

  vector<Type> artnum_t1_ll = dnorm(x_artnum_t1, A_j_t1, sd_A_j_t1, true);
  val -= sum(artnum_t1_ll);

  // Calculate model outputs

  DATA_SPARSE_MATRIX(A_out);
  DATA_SPARSE_MATRIX(A_anc_out);
  DATA_INTEGER(calc_outputs);
  DATA_INTEGER(report_likelihood)

  if(calc_outputs) {

    vector<Type> population_t1_out(A_out * population_t1);

    vector<Type> plhiv_t1_out(A_out * plhiv_t1);
    vector<Type> rho_t1_out(plhiv_t1_out / population_t1_out);

    vector<Type> artnum_t1_out(A_out * artnum_t1);
    vector<Type> alpha_t1_out(artnum_t1_out / plhiv_t1_out);
    vector<Type> artattend_t1_out(A_out * (A_artattend_mf * artnum_ij_t1));
    vector<Type> artattend_ij_t1_out(A_art_reside_attend * artnum_ij_t1);
    vector<Type> untreated_plhiv_num_t1_out(plhiv_t1_out - artnum_t1_out);

    // Calculate number of PLHIV who attend facility in district i; denominator for artattend
    vector<Type> plhiv_attend_ij_t1((Xart_idx * plhiv_t1) * (Xart_gamma * gamma_art_t1));
    vector<Type> plhiv_attend_t1_out(A_out * (A_artattend_mf * plhiv_attend_ij_t1));
    vector<Type> untreated_plhiv_attend_t1_out(plhiv_attend_t1_out - artattend_t1_out);

    vector<Type> infections_t1_out(A_out * infections_t1);
    vector<Type> lambda_t1_out(infections_t1_out / (population_t1_out - plhiv_t1_out));

    vector<Type> anc_clients_t1_out(A_anc_out * anc_clients_t1);
    vector<Type> anc_plhiv_t1_out(A_anc_out * anc_plhiv_t1);
    vector<Type> anc_already_art_t1_out(A_anc_out * anc_already_art_t1);

    // Note: assuming that:
    //  (1) anc_known_pos is equivalent to anc_already_art
    //  (2) All ANC attendees are diagnosed and initated on ART.
    vector<Type> anc_art_new_t1_out(anc_plhiv_t1_out - anc_already_art_t1_out);
    vector<Type> anc_known_pos_t1_out(anc_already_art_t1_out);
    vector<Type> anc_tested_pos_t1_out(anc_plhiv_t1_out - anc_known_pos_t1_out);
    vector<Type> anc_tested_neg_t1_out(anc_clients_t1_out - anc_plhiv_t1_out);

    vector<Type> anc_rho_t1_out(anc_plhiv_t1_out / anc_clients_t1_out);
    vector<Type> anc_alpha_t1_out(anc_already_art_t1_out / anc_plhiv_t1_out);

    REPORT(population_t1_out);
    REPORT(rho_t1_out);
    REPORT(plhiv_t1_out);
    REPORT(alpha_t1_out);
    REPORT(artnum_t1_out);
    REPORT(artattend_t1_out);
    REPORT(artattend_ij_t1_out);
    REPORT(untreated_plhiv_num_t1_out);
    REPORT(plhiv_attend_t1_out);
    REPORT(untreated_plhiv_attend_t1_out);
    REPORT(lambda_t1_out);
    REPORT(infections_t1_out);
    REPORT(anc_clients_t1_out);
    REPORT(anc_plhiv_t1_out);
    REPORT(anc_already_art_t1_out);
    REPORT(anc_art_new_t1_out);
    REPORT(anc_known_pos_t1_out);
    REPORT(anc_tested_pos_t1_out);
    REPORT(anc_tested_neg_t1_out);
    REPORT(anc_rho_t1_out);
    REPORT(anc_alpha_t1_out);

  }

  if(report_likelihood){

    REPORT(hhs_prev_ll);
    REPORT(hhs_artcov_ll);
    REPORT(artnum_t1_ll);
    REPORT(anc_rho_obs_t1_ll);

  }

  // Adam addition to include reporting of all latent field and hyper-parameter elements

  // Hyper-parameter
  REPORT(logit_phi_rho_x)
  REPORT(log_sigma_rho_x)
  REPORT(logit_phi_rho_xs)
  REPORT(log_sigma_rho_xs)
  REPORT(logit_phi_rho_a)
  REPORT(log_sigma_rho_a)
  REPORT(logit_phi_rho_as)
  REPORT(log_sigma_rho_as)
  REPORT(log_sigma_rho_xa)
  REPORT(logit_phi_alpha_x)
  REPORT(log_sigma_alpha_x)
  REPORT(logit_phi_alpha_xs)
  REPORT(log_sigma_alpha_xs)
  REPORT(logit_phi_alpha_a)
  REPORT(log_sigma_alpha_a)
  REPORT(logit_phi_alpha_as)
  REPORT(log_sigma_alpha_as)
  REPORT(log_sigma_alpha_xa)
  REPORT(OmegaT_raw)
  REPORT(log_betaT)
  REPORT(log_sigma_lambda_x)
  REPORT(log_sigma_ancrho_x)
  REPORT(log_sigma_ancalpha_x)
  REPORT(log_sigma_or_gamma)

  // Latent field
  REPORT(beta_rho)
  REPORT(beta_alpha)
  REPORT(beta_lambda)
  REPORT(beta_anc_rho)
  REPORT(beta_anc_alpha)
  REPORT(u_rho_x)
  REPORT(us_rho_x)
  REPORT(u_rho_xs)
  REPORT(us_rho_xs)
  REPORT(u_rho_a)
  REPORT(u_rho_as)
  REPORT(u_rho_xa)
  REPORT(u_alpha_x)
  REPORT(us_alpha_x)
  REPORT(u_alpha_xs)
  REPORT(us_alpha_xs)
  REPORT(u_alpha_a)
  REPORT(u_alpha_as)
  REPORT(u_alpha_xa)
  REPORT(ui_lambda_x)
  REPORT(ui_anc_rho_x)
  REPORT(ui_anc_alpha_x)
  REPORT(log_or_gamma)

  return val;
}
