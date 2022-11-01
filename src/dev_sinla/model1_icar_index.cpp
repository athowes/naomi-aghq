// model1_icar_index.cpp

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n);       // Number of regions
  DATA_INTEGER(i);       // Index i
  DATA_VECTOR(y_prev);   // Survey response
  DATA_VECTOR(m_prev);   // Survey sample size
  DATA_SPARSE_MATRIX(Q); // ICAR structure matrix
  DATA_SCALAR(Q_rank);   // Rank of Q

  // Parameter block
  PARAMETER(beta_prev);               // Survey intercept
  PARAMETER(phi_prev_i);              // Survey spatial effect for index i
  PARAMETER_VECTOR(phi_prev_minus_i); // Survey spatial effect for index -i
  PARAMETER(log_sigma_phi_prev);      // Survey log standard deviation of spatial effects

  // Transformed parameters block
  Type sigma_phi_prev = exp(log_sigma_phi_prev);

  vector<Type> phi_prev(n);
  int k = 0;
  for (int j = 0; j < n; j++) {
    if (j == i) {
      phi_prev(j) = phi_prev_i;
    } else {
      phi_prev(j) = phi_prev_minus_i(k);
      k++;
    }
  }

  vector<Type> eta_prev(beta_prev + sigma_phi_prev * phi_prev);

  // Initialise negative log-likelihood
  Type nll(0.0);

  // Priors
  nll -= dnorm(sigma_phi_prev, Type(0), Type(2.5), true) + log_sigma_phi_prev;
  nll -= dnorm(beta_prev, Type(-2), Type(1), true);

  // ICAR spatial
  nll -= Q_rank * 0.5 * log_sigma_phi_prev - 0.5 * exp(log_sigma_phi_prev) * (phi_prev * (Q * phi_prev)).sum();
  nll -= dnorm(phi_prev.sum(), Type(0.0), Type(0.001) * phi_prev.size(), true); // Soft sum-to-zero constraint

  // Likelihood
  nll -= dbinom_robust(y_prev, m_prev, eta_prev, true).sum();

  // Generated quantities block
  vector<Type> rho_prev(invlogit(eta_prev));    // Posterior prevalence estimates
  Type tau_phi_prev = 1/pow(sigma_phi_prev, 2); // Precision of spatial effects

  return(nll);
}
