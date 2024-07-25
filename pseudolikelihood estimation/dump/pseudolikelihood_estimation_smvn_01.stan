// multivariate skew-normal model with truncated-normal likelihood
functions {
  
  // Multivariate skew-normal log-density, with location vector xi,  
  // vcov matrix Omega, and slant vector alpha.
  // This function was optimized to decrease the looping computations.
  // ones is rep_vector(1.0, rows(X)) used to repeat the xi and sigma vectors 
  // in matrix form. It's passed as data to avoid creating it in every
  // iteration.
  vector multi_skew_normal_ld(data matrix X, data vector ones, row_vector xi, 
                              matrix Omega, vector alpha) {
    int N = rows(X); // X has cases in rows now
    int K = cols(X);
    vector[N] log_den;
    vector[N] log_cum;
    
    // centred and scaled x
    matrix[N, K] x_centred = X - ones * xi;
    matrix[N, K] z = x_centred ./ (ones * sqrt(diagonal(Omega))'); 
                                // just divides by sigma
    
    // Precision matrix
    matrix[K, K] P = inverse_spd(Omega);
    
    // normal density
    matrix[N, K] m_temp = x_centred * P; 
    vector[N] Q = rows_dot_product(m_temp, x_centred);
    log_den = -0.5 * Q;
    
    // cumulative probability
    log_cum = log(Phi(z * alpha)); 
    
    // skew-log-density
    return log(2) + log_den + log_cum;
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, K] X; 
  real L; // trunction for y
  
  row_vector[K] xi_lower;
  row_vector[K] xi_upper;
  real prior_xi_logit_sd;
  
  real prior_b0_mean;
  real prior_b0_sd;
  
  vector[K] prior_sigma_sd;
  real prior_sigma_nu;
  
  real prior_tau_sd;
  real prior_tau_nu;
  
  real prior_alpha_sd;
  
  real prior_corr_eta;
}

transformed data {
  vector[N] ones = rep_vector(1.0, N);
}

parameters {
  real b0_raw;
  real<lower = 0> tau; // student_t prior
  
  row_vector[K] xi_logit_raw;
  vector[K] alpha_raw;
  vector<lower = 0.01>[K] sigma; // student_t prior
 
  cholesky_factor_corr[K] Lcorr;
}

transformed parameters{
  real b0 = b0_raw * prior_b0_sd + prior_b0_mean;
  
  row_vector[K] xi_unit = inv_logit(xi_logit_raw * prior_xi_logit_sd);
  row_vector[K] xi = xi_unit .* (xi_upper - xi_lower) + xi_lower;
  
  vector[K] alpha = alpha_raw * prior_alpha_sd;
  
  corr_matrix[K] Rho; 
  cov_matrix[K] Omega;
  vector[N] mu;

  Rho = multiply_lower_tri_self_transpose(Lcorr); // = Lcorr * Lcorr'
  Omega = quad_form_diag(Rho, sigma);             // = diag_matrix(sigma) * Rho * diag_matrix(sigma)
  mu = b0 + multi_skew_normal_ld(X, ones, xi, Omega, alpha);
}

model {
  // likelihood
  for(n in 1:N)
    y[n] ~ normal(mu[n], tau)T[L, ];
  
  // priors
  b0_raw ~ std_normal();
  xi_logit_raw ~ std_normal();
  alpha_raw ~ std_normal();
  
  tau ~ student_t(prior_tau_nu, 0, prior_tau_sd);
  Lcorr ~ lkj_corr_cholesky(prior_corr_eta);
  
  for(k in 1:K)
    sigma[k] ~ student_t(prior_sigma_nu, 0, prior_sigma_sd[k]);
}