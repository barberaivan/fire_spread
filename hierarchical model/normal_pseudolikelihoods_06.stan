functions {
  
 /*
 * Evaluate the unnormalized log-density of a multivariate skew-normal variable,
 * scaled to have xi = 0, omega = 1 to simplify computations.
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
                              data matrix P, // inverse correlation matrix
                              data vector alpha) { 
    real log_den;
    real log_cum;
    
    // normal density
    real q = x' * P * x; 
    log_den = -0.5 * q; 
    
    // cumulative probability
    log_cum = log(Phi(x' * alpha));
    
    // skew-log-density
    return log_den + log_cum; // + log(2), but this is unnormalized
  }
}

data {
  int<lower=0> J;
  int<lower=0> K;
  int<lower=0> Kf;
  int<lower=0> Kr;
  int<lower=0> idfix[Kf];
  int<lower=0> idran[Kr];
  matrix[K, K] Ps[J];
  matrix[K, J] alphas;
}

parameters {
  // varying coefficients
  matrix[Kr, J] x;
  // constant coefficients
  vector[Kf] xf;
}

model {
  vector[K] xmerge; 
  for(j in 1:J) {
    xmerge[idfix] = xf;
    xmerge[idran] = x[, j];
    xmerge ~ multi_skew_normal(Ps[j], alphas[, j]);
  }
}