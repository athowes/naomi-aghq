// epil_index.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(toggle) // 0 corresponds to beta, 1 to nu, 2 to epsilon
  DATA_INTEGER(i); // Index i
  DATA_INTEGER(N);
  DATA_INTEGER(J);
  DATA_INTEGER(K);
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_MATRIX(E); // Epsilon matrix

  // Parameter block
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(nu);
  PARAMETER_VECTOR(epsilon);

  // If toggle == 0 then set-up minus i for beta
  if(toggle == 0) {

    PARAMETER(beta_i);
    PARAMETER_VECTOR(beta_minus_i);

    // vector<Type> beta(K);
    int k = 0;
    for (int j = 0; j < K; j++) {
      if (j + 1 == i) {
        beta(j) = beta_i;
      } else {
        beta(j) = beta_minus_i(k);
        k++;
      }
    }
  }

  // If toggle == 1 then set-up minus i for nu
  if(toggle == 1) {

    PARAMETER(nu_i);
    PARAMETER_VECTOR(nu_minus_i);

    // vector<Type> nu(N * J);
    int k = 0;
    for (int j = 0; j < (N * J); j++) {
      if (j + 1 == i) {
        nu(j) = nu_i;
      } else {
        nu(j) = nu_minus_i(k);
        k++;
      }
    }
  }

  // If toggle == 2 then set-up minus i for epsilon
  if(toggle == 2) {

    PARAMETER(epsilon_i);
    PARAMETER_VECTOR(epsilon_minus_i);

    // vector<Type> epsilon(N);
    int k = 0;
    for (int j = 0; j < N; j++) {
      if (j + 1 == i) {
        epsilon(j) = epsilon_i;
      } else {
        epsilon(j) = epsilon_minus_i(k);
        k++;
      }
    }
  }

  PARAMETER(l_tau_epsilon);
  PARAMETER(l_tau_nu);

  // Transformed parameters block
  Type tau_epsilon = exp(l_tau_epsilon);
  Type tau_nu = exp(l_tau_nu);
  Type sigma_epsilon = sqrt(1 / tau_epsilon);
  Type sigma_nu = sqrt(1 / tau_nu);
  vector<Type> eta(X * beta + nu + E * epsilon); // Linear predictor
  vector<Type> lambda(exp(eta));

  // Initialise negative log-likelihood
  Type nll;
  nll = Type(0.0);

  // Priors
  // Note: dgamma() is parameterised as (shape, scale); INLA parameterised as (shape, rate)
  nll -= dlgamma(l_tau_epsilon, Type(0.001), Type(1.0 / 0.001), true);
  nll -= dlgamma(l_tau_nu, Type(0.001), Type(1.0 / 0.001), true);
  nll -= dnorm(epsilon, Type(0), sigma_epsilon, true).sum();
  nll -= dnorm(nu, Type(0), sigma_nu, true).sum();
  nll -= dnorm(beta, Type(0), Type(100), true).sum();

  // Likelihood
  nll -= dpois(y, lambda, true).sum();

  // ADREPORT
  ADREPORT(tau_epsilon);
  ADREPORT(tau_nu);

  return(nll);
}
