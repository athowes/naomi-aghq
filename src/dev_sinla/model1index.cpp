// model1index.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);      // Number of regions
  DATA_VECTOR(y_prev);  // Vector of survey responses
  DATA_VECTOR(m_prev);  // Vector of survey sample sizes
  DATA_INTEGER(i);      // Index of latent field element marginal

  // Parameter block
  PARAMETER(beta_prev);              // Survey intercept
  PARAMETER(phi_prev_i);             // Survey spatial effect for index i
  PARAMETER_VECTOR(phi_prev_minus_i) // Survey spatial effect for index -i
  PARAMETER(log_sigma_phi_prev);     // Survey log standard deviation of spatial effects

  // Transformed parameters block
  Type sigma_phi_prev = exp(log_sigma_phi_prev);
  // Create phi_prev = (phi_prev_i, phi_prev_minus_i)
  vector<Type> eta_prev(beta_prev + sigma_phi_prev * phi_prev);
  // Assuming we can pick out elements of a vector in TMB then:
  // Create y_prev_i here
  // Create y_prev_minus_i here

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_prev, Type(0), Type(2.5), true) + log_sigma_phi_prev;
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);
  nll -= dnorm(phi_prev_i, Type(0), Type(1), true).sum();

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_prev, true).sum();

  // Generated quantities block
  vector<Type> rho_prev(invlogit(eta_prev));    // Posterior prevalence estimates
  Type tau_phi_prev = 1/pow(sigma_phi_prev, 2); // Precision of spatial effects

  return(nll);
}
