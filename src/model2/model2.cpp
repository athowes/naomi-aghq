// model2.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);      // Number of regions
  DATA_VECTOR(y_prev);  // Vector of survey responses
  DATA_VECTOR(m_prev);  // Vector of survey sample sizes
  DATA_VECTOR(y_anc);   // Vector of ANC response
  DATA_VECTOR(m_anc);   // Vector of ANC sample sizes

  // Parameter block
  PARAMETER(beta_prev);            // Survey intercept
  PARAMETER_VECTOR(phi_prev);      // Survey spatial effect
  PARAMETER(log_sigma_phi_prev);   // Survey log precision of spatial effects
  PARAMETER(beta_anc);             // ANC intercept
  PARAMETER_VECTOR(b_anc);         // ANC bias effect
  PARAMETER(log_sigma_b_anc);    // ANC log precision of bias effects

  // Transformed parameters block
  Type sigma_phi_prev = exp(log_sigma_phi_prev);
  vector<Type> eta_prev(beta_prev + sigma_phi_prev * phi_prev);

  Type sigma_b_anc = exp(log_sigma_b_anc);
  vector<Type> eta_anc(eta_prev + beta_anc + sigma_b_anc * b_anc);

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_prev, Type(0), Type(2.5), true) + log_sigma_phi_prev;
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);
  nll -= dnorm(phi_prev, Type(0), Type(1), true).sum();

  nll -= dnorm(sigma_b_anc, Type(0), Type(2.5), true) + log_sigma_b_anc;
  nll -= dnorm(beta_anc, Type(0), Type(1), true);
  nll -= dnorm(b_anc, Type(0), Type(1), true).sum();

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_prev, true).sum();
  nll -= dbinom_robust(y_anc, m_anc, eta_anc, true).sum();

  // Generated quantities block
  vector<Type> rho_prev(invlogit(eta_prev));
  Type tau_phi_prev = 1 / pow(sigma_phi_prev, 2);

  vector<Type> rho_anc(invlogit(eta_anc));
  Type tau_b_anc = 1 / pow(sigma_b_anc, 2);

  // ADREPORT
  ADREPORT(beta_prev);
  ADREPORT(tau_phi_prev);
  ADREPORT(phi_prev);
  ADREPORT(rho_prev);

  ADREPORT(beta_anc);
  ADREPORT(tau_b_anc);
  ADREPORT(b_anc);
  ADREPORT(rho_anc);

  return(nll);
}
