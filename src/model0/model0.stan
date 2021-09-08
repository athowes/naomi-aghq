// model0.stan

data {
  int<lower=1> n; // Number of regions
  int y_prev[n];       // Vector of responses
  int m_prev[n];       // Vector of sample sizes
}

parameters {
  real beta_rho; // Intercept
}

model {
  y_prev ~ binomial_logit(m_prev, beta_rho);
  beta_rho ~ normal(-2, 1);
}
