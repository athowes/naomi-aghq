// joint-sae.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block

  DATA_INTEGER(N_reg);
  DATA_VECTOR(pop15pl_i);
  DATA_VECTOR(pop15to49_i);
  DATA_VECTOR(art15pl_i);
  DATA_SCALAR(prev_ratio);  // Ratio of age 15+ prevalence to age 15-49 prevalence

  // Adjacency structure
  DATA_INTEGER(N_adj);
  DATA_IVECTOR(n_nb);
  DATA_IVECTOR(adj_i);
  DATA_IVECTOR(adj_j);
  DATA_IVECTOR(idx_transpose);

  // Survey prevalence
  DATA_INTEGER(Nobs_prev); // Number of observations
  DATA_IVECTOR(idx_prev);  // Region index
  DATA_VECTOR(prev_est);
  DATA_VECTOR(prev_se);

  // Survey ART coverage
  DATA_INTEGER(Nobs_arv); // Number of observations
  DATA_IVECTOR(idx_arv);  // Region index
  DATA_VECTOR(arv_est);
  DATA_VECTOR(arv_se);

  // Routine ANC testing
  DATA_INTEGER(Nobs_anc1_prev);
  DATA_IVECTOR(anc1_prev_idx);
  DATA_IVECTOR(anc1_prev_n);
  DATA_IVECTOR(anc1_prev_x);

  DATA_INTEGER(Nobs_anc1_artcov);
  DATA_IVECTOR(anc1_artcov_idx);
  DATA_IVECTOR(anc1_artcov_n);
  DATA_IVECTOR(anc1_artcov_x);

  // Flags for model components
  DATA_INTEGER(flag_prev);
  DATA_INTEGER(flag_artnum);
  DATA_INTEGER(flag_artcov);
  DATA_INTEGER(flag_ancprev);
  DATA_INTEGER(flag_ancartcov);

  // Transformed data block

  // Logit transformed survey observations
  vector<Type> l_prev_est(logit(prev_est));
  vector<Type> l_prev_se(prev_se / (prev_est * (1 - prev_est)));
  vector<Type> l_arv_est(logit(arv_est));
  vector<Type> l_arv_se(arv_se / (arv_est * (1 - arv_est)));

  vector<Type> art_attend_prior_alpha;

  // for(idx in 1:N_adj){
  //   int i;
  //   int j;
  //   i = adj_i[idx];
  //   j = adj_j[idx];
  //
  //   // alpha = 19 if same district and 0.5/[num neighbors] implies Dirichlet prior with
  //   // mean 0.94 in same region, 80% mass > 0.9
  //   if(i == j)
  //     art_attend_prior_alpha[idx] = 19;
  //   else
  //     art_attend_prior_alpha[idx] =  1.0/(n_nb[i] - 1);
  // }

  // Parameter block

  // Prevalence
  PARAMETER(l_rho0);
  PARAMETER(sigma_l_rho);
  PARAMETER_VECTOR(l_rho_i);

  PARAMETER(l_rho_ancbias0);
  PARAMETER(sigma_rho_ancbias);
  PARAMETER_VECTOR(l_rho_ancbias_raw);

  // ART coverage
  PARAMETER(l_alpha0);
  PARAMETER(sigma_l_alpha);
  PARAMETER_VECTOR(l_alpha_i_raw);

  PARAMETER(l_alpha_ancbias0);
  PARAMETER(sigma_alpha_ancbias);
  PARAMETER_VECTOR(l_alpha_ancbias_raw);

  // Distribution of district of residence for ART patients
  PARAMETER_VECTOR(pi_raw);

  // Transformed parameters block

  vector<Type> rho_i(inv_logit(l_rho_i));
  vector<Type> rho_anc_i(inv_logit(l_rho_i + l_rho_ancbias0 + sigma_rho_ancbias * l_rho_ancbias_raw));
  vector<Type> alpha_i;
  vector<Type> prop_art_i;

  vector<Type> art15pl_attend_i;
  vector<Type> pi_art;
  vector<Type> q_art;

  if (flag_artnum) {
    alpha_i = inv_logit(l_alpha0 + sigma_l_alpha * l_alpha_i_raw);
    prop_art_i = prev_ratio * rho_i .* alpha_i; // Proportion on ART
  } else {
    prop_art_i = inv_logit(l_alpha_i_raw);
    alpha_i = prop_art_i / (prev_ratio * rho_i);
  }

  int cum_nb = 0;
  // Something going on here

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors

  nll -= dnorm(l_rho0, Type(-2), Type(5), true);
  nll -= dnorm(sigma_l_rho, Type(0), Type(2.5), true);

  if (flag_prev == 1) {
    nll -= dnorm(l_rho_i, l_rho0, sigma_l_rho, true).sum();
  }

  nll -= dnorm(l_alpha0, Type(0), Type(5), true);
  nll -= dnorm(sigma_l_alpha, Type(0), Type(2.5), true);

  if (flag_artnum == 1) {
    nll -= dnorm(l_alpha_i_raw, Type(0), Type(1), true).sum();
  } else {
    nll -= -log(prop_art_i) - log1m(prop_art_i); // What's log1m?
  }

  nll -= dnorm(l_rho_ancbias0, Type(0), Type(5), true);
  nll -= dnorm(sigma_rho_ancbias, Type(0), Type(2.5), true);
  nll -= dnorm(l_rho_ancbias_raw, Type(0), Type(1), true).sum();

  nll -= dnorm(l_alpha_ancbias0, Type(0), Type(5), true);
  nll -= dnorm(sigma_alpha_ancbias, Type(0), Type(2.5), true);
  nll -= dnorm(l_alpha_ancbias_raw, Type(0), Type(1), true).sum();

  // Prior: pi_art ~ dirichlet(...);
  nll -= dlgamma(pi_raw, art_attend_prior_alpha, Type(1.0)).sum(); // Parameterisation right here?

  // Likelihood

  // Prevalence likelihood
  nll -= dnorm(l_prev_est, l_rho_i[idx_prev], l_prev_se, true).sum();

  // ANC prevalence data likelihood
  if (flag_ancprev == 1) {
    nll -= dbinom(anc1_prev_x, anc1_prev_n, rho_anc_i[anc1_prev_idx], true).sum();
  }

  // ART data likelihood
  vector<Type> pi_art_ij(prop_art_i[adj_i] * pi_art);

  // TODO

  // ADREPORT

  // TODO

  return(nll);
}
