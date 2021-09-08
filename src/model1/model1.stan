// model1.stan

data {
  int<lower=1> n; // Number of regions
  int y_prev[n];       // Vector of responses
  int m_prev[n];       // Vector of sample sizes
}

parameters {
  real beta_rho;             // Intercept
  vector[n] phi_rho;         // Spatial effects
  real<lower=0> sigma_phi_rho;   // Standard deviation of spatial effects
}

transformed parameters {
  real tau_phi_rho = 1 / sigma_phi_rho^2;                 // Precision of spatial effects
  vector[n] eta_rho = beta_rho + sigma_phi_rho * phi_rho; // Linear predictor
}

model {
  y_prev ~ binomial_logit(m_prev, eta_rho);
  phi_rho ~ normal(0, 1); // phi has variance one
  beta_rho ~ normal(-2, 1);
  sigma_phi_rho ~ normal(0, 2.5);
}
