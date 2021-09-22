// model1.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);      // Number of regions
  DATA_VECTOR(y_prev);  // Vector of survey responses
  DATA_VECTOR(m_prev);  // Vector of survey sample sizes

  // Parameter block
  PARAMETER(beta_prev);            // Survey intercept
  PARAMETER_VECTOR(phi_prev);      // Survey spatial effect
  PARAMETER(log_sigma_phi_prev);   // Survey log standard deviation of spatial effects

  // Transformed parameters block
  Type sigma_phi_prev = exp(log_sigma_phi_prev);
  vector<Type> eta_prev(beta_prev + sigma_phi_prev * phi_prev);

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_prev, Type(0), Type(2.5), true) + log_sigma_phi_prev;
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);
  nll -= dnorm(phi_prev, Type(0), Type(1), true).sum();

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_prev, true).sum();

  // Generated quantities block
  vector<Type> rho_prev(invlogit(eta_prev));    // Posterior prevalence estimates
  Type tau_phi_prev = 1/pow(sigma_phi_prev, 2); // Precision of spatial effects

  // ADREPORT
  ADREPORT(beta_prev);
  ADREPORT(tau_phi_prev);
  ADREPORT(phi_prev);
  ADREPORT(rho_prev);

  return(nll);
}
