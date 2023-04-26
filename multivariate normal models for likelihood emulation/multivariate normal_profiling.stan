// multivariate normal model

/*
This code tries many ways to compute the unnormalized multivariate normal
density. 
The 10th function is the best (same computation as that in the sn package)


USAR ESTOOO
mdivide_right_spd 

*/


functions {
  /*
  // Multivariate skew-normal log-density, with location vector xi, and 
  // vcov matrix Omega
  vector multi_normal_ld(data matrix X, vector xi, matrix Omega) {
    int N = cols(X); // X has cases in columns, to avoid transposing and looping
                     // over columns
    int K = rows(X);
    vector[N] log_den;
    vector[K] x_centred;
    matrix[K, K] P = inverse_spd(Omega);

    for(n in 1:N) {
      x_centred = X[, n] - xi;
      log_den[n] = -0.5 * x_centred' * P * x_centred;
    }
    
    return log_den;
  }
  */
  
  // A few ways to compute the loop inside multi_normal_ld (the heavy step)
  
  vector ld_01(data matrix X, vector xi, matrix P) {
    int N = cols(X); // X has cases in columns, to avoid transposing and looping
                     // over columns
    int K = rows(X);
    vector[N] log_den;
    vector[K] x_centred;

    for(n in 1:N) {
      x_centred = X[, n] - xi;
      log_den[n] = -0.5 * x_centred' * P * x_centred;
    }
    
    return log_den;
  }
  
  vector ld_02(data matrix X_long, vector xi, matrix P) {
    int N = rows(X_long); // X has cases in rows now
    int K = cols(X_long);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    vector[N] ones = rep_vector(1, N);
    matrix[N, K] Xi = ones * xi';
    matrix[N, K] x_centred = X_long - Xi;
    
    // firt multiplication
    matrix[N, K] m1 = x_centred * P;
    // second (giant) multiplication
    matrix[N, N] m2 = m1 * x_centred';
    
    // result
    log_den = -0.5 * diagonal(m2);
    return log_den;
  }
  
  vector ld_03(data matrix X_long, vector xi, matrix P) {
    int N = rows(X_long); // X has cases in rows now
    int K = cols(X_long);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    vector[N] ones = rep_vector(1, N);
    matrix[N, K] Xi = ones * xi';
    matrix[N, K] x_centred = X_long - Xi;
    
    // firt multiplication
    matrix[N, K] m1 = x_centred * P;
    // second (giant) multiplication
    vector[N] temp;
    for(n in 1:N) {
      temp[n] = m1[n, ] * x_centred[n, ]';
    }
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_04(data matrix X_long, vector xi, matrix P) {
    int N = rows(X_long); // X has cases in rows now
    int K = cols(X_long);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    vector[N] ones = rep_vector(1, N);
    matrix[N, K] Xi = ones * xi';
    matrix[N, K] x_centred = X_long - Xi;
    
    // firt multiplication
    matrix[N, K] m1 = x_centred * P;
    // second (giant) multiplication
    vector[N] temp = rows_dot_product(m1, x_centred);
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }

  vector ld_05(data matrix X_long, vector xi, matrix P) {
    int N = rows(X_long); // X has cases in rows now
    int K = cols(X_long);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    vector[N] ones = rep_vector(1, N);
    matrix[N, K] Xi = ones * xi';
    matrix[N, K] x_centred = X_long - Xi;
    matrix[K, N] x_centred_t = x_centred';
    
    // firt multiplication
    matrix[K, N] m1 = (x_centred * P)';
    // second (giant) multiplication
    vector[N] temp;
    for(n in 1:N) {
      temp[n] = dot_product(m1[, n], x_centred_t[, n]);
    }
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_06(data matrix X_long, vector xi, matrix P) {
    int N = rows(X_long); // X has cases in rows now
    int K = cols(X_long);
    vector[N] log_den;
    // matrix[K, N] x_centred;
    vector[N] ones = rep_vector(1, N);
    matrix[N, K] Xi = ones * xi';
    matrix[N, K] x_centred = X_long - Xi;
    matrix[K, N] x_centred_t = x_centred';
    
    // firt multiplication
    matrix[K, N] m1 = P * x_centred_t;
    // second (giant) multiplication
    vector[N] temp;
    for(n in 1:N) {
      temp[n] = dot_product(m1[, n], x_centred_t[, n]);
    }
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_07(data matrix X, data row_vector ones, vector xi, matrix P) {
    int K = rows(X); // X has cases in columns
    int N = cols(X);
    vector[N] log_den;
    
    // get x_centred matrix
    matrix[K, N] Xi = xi * ones;
    matrix[K, N] x_centred = X - Xi;

    // firt multiplication
    matrix[K, N] m1 = P * x_centred;
    // second (giant) multiplication
    vector[N] temp;
    for(n in 1:N) {
      temp[n] = dot_product(m1[, n], x_centred[, n]);
    }
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_08(data matrix X, data row_vector ones, vector xi, matrix P) {
    int K = rows(X); // X has cases in columns
    int N = cols(X);
    vector[N] log_den;
    
    // get x_centred matrix
    matrix[K, N] Xi = xi * ones;
    matrix[K, N] x_centred = X - Xi;

    // first multiplication
    matrix[K, N] m_temp = (P * x_centred) .* x_centred; // note that P = P'
    // second (giant) multiplication
    vector[N] temp;
    for(n in 1:N) {
      temp[n] = sum(m_temp[, n]);
    }
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_09(data matrix X, data row_vector ones_N, data row_vector ones_K, 
               vector xi, matrix P) {
    int K = rows(X); // X has cases in columns
    int N = cols(X);
    vector[N] log_den;
    
    // get x_centred matrix
    matrix[K, N] Xi = xi * ones_N;
    matrix[K, N] x_centred = X - Xi;

    // first multiplication
    matrix[K, N] m_temp = (P * x_centred) .* x_centred; // note that P = P'
    // second (giant) multiplication
    vector[N] temp = (ones_K * m_temp)';
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_10(data matrix X_long, data vector ones_N_vec, data vector ones_K_vec, 
               row_vector xi, matrix Omega) {
    
    int K = cols(X_long); // X has cases in columns
    int N = rows(X_long);
    vector[N] log_den;
    matrix[K, K] P = inverse_spd(Omega);
    
    // get x_centred matrix
    matrix[N, K] Xi = ones_N_vec * xi;
    matrix[N, K] x_centred = X_long - Xi;

    // first multiplication
    matrix[N, K] m_temp = (x_centred * P) .* x_centred; // note that P = P'
    vector[N] temp = m_temp * ones_K_vec;
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_11(data matrix X_long, data vector ones_N_vec, data vector ones_K_vec, 
               row_vector xi, matrix Omega) {
    
    int K = cols(X_long); // X has cases in columns
    int N = rows(X_long);
    vector[N] log_den;

    // get x_centred matrix
    matrix[N, K] Xi = ones_N_vec * xi;
    matrix[N, K] x_centred = X_long - Xi;

    // first multiplication
    matrix[N, K] m_temp = mdivide_right_spd(x_centred, Omega) .* x_centred; 
    vector[N] temp = m_temp * ones_K_vec;
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  // I think this is how Simon Wood computes it, but the results do not match
  vector ld_12(data matrix X_long, data vector ones_N_vec, data vector ones_K_vec, 
               row_vector xi, matrix Omega) {
    
    int K = cols(X_long); // X has cases in columns
    int N = rows(X_long);
    vector[N] log_den;

    // get x_centred matrix
    matrix[N, K] Xi = ones_N_vec * xi;
    matrix[N, K] x_centred = X_long - Xi;
    
    // Choleski factor for Omega
    matrix[K, K] L = cholesky_decompose(Omega);

    // first multiplication
    matrix[N, K] m_temp = mdivide_right_tri_low(x_centred, Omega); 
    vector[N] temp = (m_temp .* m_temp) * ones_K_vec;
    
    // result
    log_den = -0.5 * temp;
    return log_den;
  }
  
  vector ld_13(data matrix X_long, data vector ones_N_vec, data vector ones_K_vec, 
               row_vector xi, matrix Omega) {
    
    int K = cols(X_long); // X has cases in columns
    int N = rows(X_long);
    vector[N] log_den;

    // get x_centred matrix
    matrix[N, K] Xi = ones_N_vec * xi;
    matrix[N, K] x_centred = X_long - Xi;

    // first multiplication
    matrix[N, K] m_temp = (x_centred / Omega) .* x_centred; 
    vector[N] temp = m_temp * ones_K_vec;
    
    // result
    log_den = -0.5 * temp;
    return log_den;
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
}

transformed data {
  matrix[N, K] X_long = X';
  row_vector[N] ones = rep_row_vector(1, N);
  row_vector[N] ones_N = rep_row_vector(1, N);
  row_vector[K] ones_K = rep_row_vector(1, K);
  vector[N] ones_N_vec = rep_vector(1, N);
  vector[K] ones_K_vec = rep_vector(1, K);
}

parameters {
  vector[K] xi;
  vector<lower = 0.05>[K] sigma; // lower > 0 to ensure numerical positive-definiteness
  cholesky_factor_corr[K] Lcorr;
  real<lower = 0> tau;
}

transformed parameters{
  corr_matrix[K] Rho; 
  cov_matrix[K] Omega;
  vector[N] mu;
  vector[N] mu_02;
  vector[N] mu_03;
  vector[N] mu_04;
  vector[N] mu_05;
  vector[N] mu_06;
  vector[N] mu_07;
  vector[N] mu_08;
  vector[N] mu_09;
  vector[N] mu_10;
  vector[N] mu_11;
  vector[N] mu_12;
  vector[N] mu_13;
  
  // definitions for mu_ld
  vector[N] log_den;
  vector[K] x_centred;
  matrix[K, K] P;
  
  Rho = multiply_lower_tri_self_transpose(Lcorr); // = Lcorr * Lcorr'
  Omega = quad_form_diag(Rho, sigma);             // = diag_matrix(sigma) * Rho * diag_matrix(sigma)
  
  // mu = multi_normal_ld(X, xi, Omega);
  
  // COMPUTE MU BY HAND HERE to profile its inner parts
  
  P = inverse_spd(Omega);

  profile("ld_01") {
    mu = ld_01(X, xi, P);
  }
  
  profile("ld_02") {
    mu_02 = ld_02(X_long, xi, P);
  }
  
  profile("ld_03") {
    mu_03 = ld_03(X_long, xi, P);
  }
  
  profile("ld_04") {
    mu_04 = ld_04(X_long, xi, P);
  }
  
  profile("ld_05") {
    mu_05 = ld_05(X_long, xi, P);
  }
  
  profile("ld_06") {
    mu_06 = ld_06(X_long, xi, P);
  }
  
  profile("ld_07") {
    mu_07 = ld_07(X, ones, xi, P);
  }

  profile("ld_08") {
    mu_08 = ld_08(X, ones, xi, P);
  }
  
  profile("ld_09") {
    mu_09 = ld_09(X, ones_N, ones_K, xi, P);
  }
  
  profile("ld_10") {
    mu_10 = ld_10(X_long, ones_N_vec, ones_K_vec, xi', Omega);
  }
  
  profile("ld_11") {
    mu_11 = ld_11(X_long, ones_N_vec, ones_K_vec, xi', Omega);
  }
  
  profile("ld_12") {
    mu_12 = ld_12(X_long, ones_N_vec, ones_K_vec, xi', Omega);
  }
  
  profile("ld_13") {
    mu_13 = ld_13(X_long, ones_N_vec, ones_K_vec, xi', Omega);
  }
}

model {
  // likelihood
  y ~ normal(mu, tau);
  
  // priors
  xi ~ normal(0, xi_prior_sd);
  sigma ~ student_t(sigma_prior_nu, 0, sigma_prior_sd); // nu, mu, sigma
  tau ~ student_t(tau_prior_nu, 0, tau_prior_sd);
  Lcorr ~ lkj_corr_cholesky(1);
}