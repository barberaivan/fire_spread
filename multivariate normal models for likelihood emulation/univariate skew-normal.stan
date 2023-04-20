// univariate skew-normal model
functions{
  /*
   * univariate skew-normal unnormalized log-density 
   */ 
  vector skew_normal_ld(data vector x, real xi, real sigma, real alpha) {
    
    int N = rows(x);
    vector[N] x_centred = x - xi;
    vector[N] shifted = alpha * x_centred / sigma;
    real prec = 1 / sigma ^ 2;
    
    vector[N] log_den = -0.5 * prec * x_centred .* x_centred; // normal_lupdf doesn't work
    vector[N] log_cum = log(0.5) + log(1 + erf(shifted / sqrt(2))); // std_normal_lcdf doesnt work
    
    return log(2) + log_den + log_cum;
  }

}

data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real xi_prior_sd;
  real sigma_prior_sd;
  real sigma_prior_nu;
  real tau_prior_sd;
  real tau_prior_nu;
  real alpha_prior_sd;
}

parameters {
  real xi;
  real<lower = 0> sigma;
  real<lower = 0> tau;
  real alpha;
}

transformed parameters {
  vector[N] mu = skew_normal_ld(x, xi, sigma, alpha);
}

model {
  // likelihood
  y ~ normal(mu, tau);
  
  // priors
  xi ~ normal(0, xi_prior_sd);
  tau ~ student_t(tau_prior_nu, 0, tau_prior_sd);
  
  alpha ~ normal(0, alpha_prior_sd);
  //sigma ~ student_t(sigma_prior_nu, 0, sigma_prior_sd); // nu, mu, sigma
  sigma ~ normal(0, sigma_prior_sd);
}