functions {
  
  // Original function, designed to evaluate the density of a data matrix
  // X:
  
  // Multivariate skew-normal log-density, with location vector xi,  
  // vcov matrix Omega, and slant vector alpha.
  // This function was optimized to decrease the looping computations.
  // ones is rep_vector(1.0, rows(X)) used to repeat the xi and sigma vectors 
  // in matrix form. It's passed as data to avoid creating it in every
  // iteration.
  /*
  vector multi_skew_normal_lpdf(matrix X, data vector ones, data row_vector xi, 
                                data matrix Omega, data vector alpha) {
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
  */

//
  
 /*
 * Evaluate the unnormalized log-density of a multivariate skew-normal variable.
 * It is designed to evaluate the density at a given vector, not matrix. 
 * As in this application, the parameters will be treated as known, they are
 * declared as data, and instead of providing Omega, their transformed 
 * quantities are (P, sigma_inverse).
 
 * @param vector x: where to evaluate log-density.
 * @param vector xi: location (~mu).
 * @param matrix P: inv(Omega).
 * @param vector alpha: slant parameter.
 * @param vector sigma_inverse: 1 / sqrt(diag(Omega)).
 
 * @return log-density at x.
 */
  real multi_skew_normal_lpdf(vector x, 
                              data vector xi, 
                              data matrix P, 
                              data vector alpha,
                              data vector sigma_inverse) { 
    int K = size(xi);
    real log_den;
    real log_cum;
    
    // centred and scaled x
    vector[K] x_centred = x - xi;
    vector[K] z = x_centred .* sigma_inverse;
    
    // normal density
    real q = x_centred' * P * x_centred; 
    log_den = -0.5 * q; 
    
    // cumulative probability
    log_cum = log(Phi(z' * alpha));
    
    // skew-log-density
    return log_den + log_cum; // + log(2), but this is unnormalized
  }
}

data {
  int<lower=0> J;
  int<lower=0> K;
  matrix[K, J] xis;
  matrix[K, K] Ps[J];
  matrix[K, J] alphas;
  matrix[K, J] sigmas_inv;
}

parameters {
  matrix[K, J] x;
}

model {
  for(j in 1:J)
    x[, j] ~ multi_skew_normal(xis[, j], Ps[j], alphas[, j], sigmas_inv[, j]);
}