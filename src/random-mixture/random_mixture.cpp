// random_mixture.cpp
// From https://github.com/kaskr/laplace_accuracy/blob/master/random_mixture/random_mixture.cpp

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(N);

  // Latent variables
  PARAMETER_MATRIX(u);   // Rows are observations, columns are groups

  // Parameters
  PARAMETER_VECTOR(mu);     // Mean by group
  PARAMETER_VECTOR(logsd);  // Sd by group

  Type ans = 0;

  // u(i,j) ~ N(mu(j), sd(j))
  vector<Type> sd = exp(logsd);
  for (int i=0; i<u.rows(); i++) {
    for (int j=0; j<u.cols(); j++) {
      ans -= dnorm(u(i, j), mu(j), sd(j), true);
      SIMULATE{ u(i, j) = rnorm(mu(j), sd(j)); }
    }
  }

  // Intensity by observation
  vector<Type> lambda = u.array().exp().rowwise().sum();
  ans -= dpois(N, lambda, true).sum();
  SIMULATE{ N = rpois(lambda); REPORT(u); REPORT(N); }

  return ans;
}
