// multivariate normal model with truncated-normal likelihood
functions {
  
  // Multivariate skew-normal log-density, with location vector xi, and 
  // vcov matrix Omega.
  // This function was optimized to decrease the looping computations.
  vector multi_normal_ld(data matrix X, data vector ones, row_vector xi, 
                         matrix Omega) {
    int N = rows(X); // X has cases in rows now
    int K = cols(X);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    matrix[N, K] Xi = ones * xi;
    matrix[N, K] x_centred = X - Xi;
    
    // Precision matrix
    matrix[K, K] P = inverse_spd(Omega);
    
    // firt multiplication
    matrix[N, K] m_temp = x_centred * P;
    // second multiplication
    vector[N] Q = rows_dot_product(m_temp, x_centred);
    
    // result
    log_den = -0.5 * Q;
    return log_den;
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[N, K] X; 
  
  real xi_prior_sd;
  real sigma_prior_sd;
  real sigma_prior_nu;
  real tau_prior_sd;
  real tau_prior_nu;
  
  real sigma_lower;
  real y_lower;
}

transformed data {
  vector[N] ones = rep_vector(1.0, N);
}

parameters {
  row_vector[K] xi;
  vector<lower = sigma_lower>[K] sigma;
  cholesky_factor_corr[K] Lcorr;
  real<lower = 0> tau;
}

transformed parameters{
  corr_matrix[K] Rho; 
  cov_matrix[K] Omega;
  vector[N] mu;

  Rho = multiply_lower_tri_self_transpose(Lcorr); // = Lcorr * Lcorr'
  Omega = quad_form_diag(Rho, sigma);             // = diag_matrix(sigma) * Rho * diag_matrix(sigma)
  mu = multi_normal_ld(X, ones, xi, Omega);
}

model {
  // likelihood
  for(n in 1:N) 
    y[n] ~ normal(mu[n], tau)T[y_lower, ];
  
  // priors
  xi ~ normal(0, xi_prior_sd);
  sigma ~ student_t(sigma_prior_nu, 0, sigma_prior_sd); // nu, mu, sigma
  tau ~ student_t(tau_prior_nu, 0, tau_prior_sd);
  Lcorr ~ lkj_corr_cholesky(1);
}