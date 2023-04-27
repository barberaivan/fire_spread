// univariate skew-normal model
functions{
  /*
   * univariate skew-normal unnormalized log-density. Important, use Phi() to 
   * compute the standard normal cdf. Computing it by hand caused problems.
   */ 
  vector skew_normal_ld(data vector x, real xi, real sigma, real alpha) {
    
    int N = rows(x);
    vector[N] x_centred = x - xi;
    vector[N] shifted = alpha * x_centred / sigma;
    real prec = 1 / sigma ^ 2;
    
    vector[N] log_den = -0.5 * prec * x_centred .* x_centred; 
    vector[N] log_cum = normal_lcdf(shifted | 0, 1); //log(Phi(shifted));
                      // change not tested
    
    return log(2) + log_den + log_cum;
  }
  
  // return all intermediate steps (for testing mu computation in R and Stan)
  matrix skew_normal_ld_debug(data vector x, real xi, real sigma, real alpha) {
    
    int N = rows(x);
    vector[N] x_centred = x - xi;
    vector[N] shifted = alpha * x_centred / sigma;
    real prec = 1 / sigma ^ 2;
    
    vector[N] log_den = -0.5 * prec * x_centred .* x_centred; // normal_lupdf doesn't work
    vector[N] log_cum = log(0.5) + log(1 + erf(shifted / sqrt(2))); // std_normal_lcdf doesnt work
    
    matrix[N, 10] result;
    result[, 1] = log_den;
    result[, 2] = log_cum;
    result[, 3] = log(2) + log_den + log_cum;
    result[, 4] = x_centred;
    result[, 5] = x_centred .* x_centred;
    result[, 6] = rep_vector(prec, N);
    result[, 7] = shifted;
    result[, 8] = erf(shifted / sqrt(2));
    result[, 9] = Phi(shifted);
    result[, 10] = log(2) + log_den + log(Phi(shifted));
    return result;
  }

}

data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  
  // priors
  real xi_prior_sd;
  real sigma_prior_sd;
  real sigma_prior_nu;
  real tau_prior_sd;
  real tau_prior_nu;
  real alpha_prior_sd;
  
  // true parameters
  real xi_true;
  real sigma_true;
  real alpha_true;
  real tau_true;
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

generated quantities {
  // true mu
  matrix[N, 10] mu_true_etc = skew_normal_ld_debug(x, xi_true, sigma_true, alpha_true);
}