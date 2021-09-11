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
  PARAMETER(beta_prev); // Intercept

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, beta_prev, true).sum();

  return(nll);
}
