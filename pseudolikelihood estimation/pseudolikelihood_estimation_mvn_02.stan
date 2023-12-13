// multivariate normal model with truncated-normal likelihood
functions {

  // Multivariate normal log-density, with location vector xi, and
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

  real prior_corr_eta;
}

transformed data {
  vector[N] ones = rep_vector(1.0, N);
}

parameters {
  real b0_raw;
  real<lower = 0> tau; // student_t prior

  row_vector[K] xi_logit_raw;
  vector<lower = 0.01>[K] sigma; // student_t prior

  cholesky_factor_corr[K] Lcorr;
}

transformed parameters{
  real b0 = b0_raw * prior_b0_sd + prior_b0_mean;

  row_vector[K] xi_unit = inv_logit(xi_logit_raw * prior_xi_logit_sd);
  row_vector[K] xi = xi_unit .* (xi_upper - xi_lower) + xi_lower;

  corr_matrix[K] Rho;
  cov_matrix[K] Omega;
  vector[N] mu;

  Rho = multiply_lower_tri_self_transpose(Lcorr); // = Lcorr * Lcorr'
  Omega = quad_form_diag(Rho, sigma);             // = diag_matrix(sigma) * Rho * diag_matrix(sigma)
  mu = b0 + multi_normal_ld(X, ones, xi, Omega);
}

model {
  // likelihood
  y ~ normal(mu, tau);

  // priors
  b0_raw ~ std_normal();
  xi_logit_raw ~ std_normal();

  tau ~ student_t(prior_tau_nu, 0, prior_tau_sd);
  Lcorr ~ lkj_corr_cholesky(prior_corr_eta);

  for(k in 1:K)
    sigma[k] ~ student_t(prior_sigma_nu, 0, prior_sigma_sd[k]);
}