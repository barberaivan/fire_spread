// multivariate skew-normal model
functions {
  
  // Multivariate skew-normal log-density, with location vector xi, vcov matrix 
  // Omega and slant vector alpha
  vector multi_skew_normal_ld(data matrix X, vector xi, matrix Omega, 
                              row_vector alpha) {
    int N = cols(X); // X has cases in columns, to avoid transposing and looping
                     // over columns
    int K = rows(X);
    vector[N] log_den;
    vector[N] log_cum;
    vector[K] x_centred;
    vector[K] sigma = sqrt(diagonal(Omega)); 
    matrix[K, K] P = inverse_spd(Omega);
                                     
    for(n in 1:N) {
      x_centred = X[, n] - xi;
      
      // normal log_density
      log_den[n] = -0.5 * x_centred' * P * x_centred;
      
      // normal log_cumulative_prob for shifted value
      log_cum[n] = log(Phi(alpha * (x_centred ./ sigma)));
    }
    
    return log(2) + log_den + log_cum;
  }
}

data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  matrix[K, N] X; // transposed to loop over columns, which is more efficient
  
  real xi_prior_sd;
  real sigma_prior_sd;
  real sigma_prior_nu;
  real tau_prior_sd;
  real tau_prior_nu;
  real alpha_prior_sd;
}

parameters {
  vector[K] xi;
  vector<lower = 0>[K] sigma;
  row_vector[K] alpha;
  cholesky_factor_corr[K] Lcorr;
  real<lower = 0> tau;
}

transformed parameters{
  corr_matrix[K] Rho; 
  cov_matrix[K] Omega;
  cov_matrix[K] P;     // precision matrix
  vector[N] mu;

  Rho = multiply_lower_tri_self_transpose(Lcorr); // = Lcorr * Lcorr'
  Omega = quad_form_diag(Rho, sigma);             // = diag_matrix(sigma) * Rho * diag_matrix(sigma)
  P = inverse_spd(Omega);
  mu = multi_skew_normal_ld(X, xi, P, alpha);
}

model {
  # likelihood
  y ~ normal(mu, tau);
  
  # priors
  xi ~ normal(0, xi_prior_sd);
  sigma ~ student_t(sigma_prior_nu, 0, sigma_prior_sd); // nu, mu, sigma
  alpha ~ normal(0, alpha_prior_sd);
  tau ~ student_t(tau_prior_nu, 0, tau_prior_sd);
  Lcorr ~ lkj_corr_cholesky(1);
}