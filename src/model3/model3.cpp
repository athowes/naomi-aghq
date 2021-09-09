// model3.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);      // Number of regions
  DATA_VECTOR(y_prev);  // Vector of survey responses
  DATA_VECTOR(m_prev);  // Vector of survey sample sizes
  DATA_VECTOR(A_art);   // Number of people on ART
  DATA_VECTOR(N_art);   // Total population size

  // Parameter block
  PARAMETER(beta_prev);            // Survey intercept
  PARAMETER_VECTOR(phi_prev);      // Survey spatial effect
  PARAMETER(log_sigma_phi_prev);   // Survey log standard deviation of spatial effects
  PARAMETER(beta_art);             // ART intercept
  PARAMETER_VECTOR(phi_art);       // ART random effects
  PARAMETER(log_sigma_phi_art);     // ART log standard deviation of random effects

  // Transformed parameters block
  Type sigma_phi_prev = exp(log_sigma_phi_prev);
  vector<Type> eta_prev(beta_prev + sigma_phi_prev * phi_prev);

  Type sigma_phi_art = exp(log_sigma_phi_art);
  vector<Type> eta_art(beta_art + sigma_phi_art * phi_art);

  vector<Type> rho_prev(invlogit(eta_prev));
  vector<Type> alpha_art(invlogit(eta_art));

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_prev, Type(0), Type(2.5), true) + log_sigma_phi_prev;
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);
  nll -= dnorm(phi_prev, Type(0), Type(1), true).sum();

  nll -= dnorm(sigma_phi_art, Type(0), Type(2.5), true) + log_sigma_phi_art;
  nll -= dnorm(beta_art, Type(0), Type(1), true);
  nll -= dnorm(phi_art, Type(0), Type(1), true).sum();

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_prev, true).sum();
  log(alpha_art) = log(A_art) + log(rho_prev) - log(N_art); // alpha = A * rho / N

  // Generated quantities block
  Type tau_phi_prev = 1/pow(sigma_phi_prev, 2);
  Type tau_phi_art = 1/pow(sigma_phi_art, 2);

  // ADREPORT
  ADREPORT(beta_prev);
  ADREPORT(tau_phi_prev);
  ADREPORT(phi_prev);
  ADREPORT(rho_prev);

  ADREPORT(beta_art);
  ADREPORT(tau_phi_art);
  ADREPORT(phi_art);
  ADREPORT(alpha_art);

  return(nll);
}
