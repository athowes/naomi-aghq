// blangiardo.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_VECTOR(y);

  // Parameter block
  PARAMETER(x);
  PARAMETER(l_theta);

  // Transformed parameters block
  Type sigma = sqrt(1 / exp(l_theta));

  Type nll;
  nll = Type(0.0);

  // Priors
  nll -= dlgamma(l_theta, Type(1.6), Type(1 / 0.4), true);
  nll -= dnorm(x, Type(-3.0), Type(2.0), true);

  // Likelihood
  nll -= dnorm(y, x, sigma, true).sum();

  // ADREPORT
  ADREPORT(x)
  ADREPORT(sigma);

  return(nll);
}
