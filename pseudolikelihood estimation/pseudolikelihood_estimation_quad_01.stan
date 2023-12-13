// multivariate normal model with truncated-normal likelihood.
// This model assumes no correlation at all.

functions {
  // Multivariate normal log-density, with location vector xi, and
  // precision matrix P.
  // This function was optimized to decrease the looping computations.
  vector multi_normal_ld(data matrix X, data vector ones, row_vector xi,
                         matrix P) {
    int N = rows(X); // X has cases in rows now
    int K = cols(X);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    matrix[N, K] Xi = ones * xi;
    matrix[N, K] x_centred = X - Xi;

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
  int truncated; // if 1, use truncated likelihood
  real L; // trunction for y

  row_vector[K] xi_lower;
  row_vector[K] xi_upper;
  real prior_xi_logit_sd;

  real prior_b0_mean;
  real prior_b0_sd;

  vector[K] prior_prec_sd;
  real prior_prec_nu;

  real prior_tau_sd;
  real prior_tau_nu;
}

transformed data {
  vector[N] ones = rep_vector(1.0, N);
}

parameters {
  real b0_raw;
  real<lower = 0> tau; // student_t prior

  row_vector[K] xi_logit_raw;
  vector<lower = 0>[K] precision; // student_t prior
}

transformed parameters{
  real b0 = b0_raw * prior_b0_sd + prior_b0_mean;

  row_vector[K] xi_unit = inv_logit(xi_logit_raw * prior_xi_logit_sd);
  row_vector[K] xi = xi_unit .* (xi_upper - xi_lower) + xi_lower;

  matrix[K, K] P = diag_matrix(precision);
  vector[N] mu;

  mu = b0 + multi_normal_ld(X, ones, xi, P);
}

model {
  // likelihood
  if(truncated == 1) {
    for(n in 1:N)
    y[n] ~ normal(mu[n], tau)T[L, ];
  } else {
    y ~ normal(mu, tau);
  }

  // priors
  b0_raw ~ std_normal();
  xi_logit_raw ~ std_normal();

  tau ~ student_t(prior_tau_nu, 0, prior_tau_sd);

  for(k in 1:K)
    precision[k] ~ student_t(prior_prec_nu, 0, prior_prec_sd[k]);
}