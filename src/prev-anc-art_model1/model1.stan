// model1.stan

data {
  int<lower=1> n;      // Number of regions
  int y_prev[n];       // Vector of responses
  int m_prev[n];       // Vector of sample sizes
}

parameters {
  real beta_prev;                 // Intercept
  vector[n] phi_prev;             // Spatial effects
  real<lower=0> sigma_phi_prev;   // Standard deviation of spatial effects
}

transformed parameters {
  real tau_phi_prev = 1 / sigma_phi_prev^2;                   // Precision of spatial effects
  vector[n] eta_prev = beta_prev + sigma_phi_prev * phi_prev; // Linear predictor
}

model {
  y_prev ~ binomial_logit(m_prev, eta_prev);
  phi_prev ~ normal(0, 1); // phi has variance one
  beta_prev ~ normal(-2, 1);
  sigma_phi_prev ~ normal(0, 2.5);
}
