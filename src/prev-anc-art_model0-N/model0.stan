// model0.stan

data {
  int<lower=1> n; // Number of regions
  int y_prev[n];  // Vector of responses
  int m_prev[n];  // Vector of sample sizes
}

parameters {
  real beta_prev; // Intercept
}

model {
  y_prev ~ binomial_logit(m_prev, beta_prev);
  beta_prev ~ normal(-2, 1);
}
