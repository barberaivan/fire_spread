// univariate normal model
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real xi_prior_sd;
  real sigma_prior_sd;
  real sigma_prior_nu;
  real tau_prior_sd;
  real tau_prior_nu;
}

parameters {
  real xi;
  real<lower = 0> sigma;
  real<lower = 0> tau;
}

transformed parameters{
  real<lower = 0> prec = 1 / sigma ^ 2;
}

model {
  vector[N] mu = -0.5 * prec * (x - xi) .* (x - xi); // does not allow a ^2 on a vector
  y ~ normal(mu, tau);
  
  # priors
  xi ~ normal(0, xi_prior_sd);
  sigma ~ student_t(sigma_prior_nu, 0, sigma_prior_sd); // nu, mu, sigma
  tau ~ student_t(tau_prior_nu, 0, tau_prior_sd);
}