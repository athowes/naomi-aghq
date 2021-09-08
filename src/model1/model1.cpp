// model1.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);      // Number of regions
  DATA_VECTOR(y_prev);  // Vector of responses
  DATA_VECTOR(m_prev);  // Vector of sample sizes

  // Parameter block
  PARAMETER(beta_rho);            // Intercept
  PARAMETER_VECTOR(phi_rho);      // Spatial effect
  PARAMETER(log_sigma_phi_rho);   // Log precision of spatial effects

  // Transformed parameters block
  Type sigma_phi_rho = exp(log_sigma_phi_rho); // Standard deviation of spatial effects
  vector<Type> eta_rho(beta_rho + sigma_phi_rho * phi_rho);

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_rho, Type(0), Type(2.5), true) + log_sigma_phi_rho;
  nll -= dnorm(beta_rho, Type(-2), Type(1), true);
  nll -= dnorm(phi_rho, Type(0), Type(1), true).sum();

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_rho, true).sum();

  // Generated quantities block
  vector<Type> rho(invlogit(eta_rho)); // Posterior prevalence estimates
  Type tau_phi_rho = 1/pow(sigma_phi_rho, 2); // Precision of spatial effects

  // ADREPORT
  ADREPORT(beta_rho);
  ADREPORT(tau_phi_rho);
  ADREPORT(phi_rho);
  ADREPORT(rho);

  return(nll);
}
