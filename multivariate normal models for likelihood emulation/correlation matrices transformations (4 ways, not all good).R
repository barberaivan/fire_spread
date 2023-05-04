build_corr <- function(y, d) {
  # test
  d <- 4
  y <- rnorm(cor_num(d))
  ##
  # upper triangular matrix of Canonical Partial Correlations
  z <- matrix(0, d, d)
  z[upper.tri(z)] <- tanh(y)
  
  # compute omega in many ways to test:
  
  # 0 and 1 are the Stan's efficient way, but they don't result in a 1-diagonal
  # correlation matrix. 2 and 3 are the fool way, but they work.
  
  ### 0
  w0 <- matrix(0, d, d)
  system.time(
    {
      w0[1, ] <- c(1, z[1, 2:d])
      
      for(i in 2:d) {
        for(j in i:d) {
          w0[i, j] <- z[i, j] * w0[i-1, j] * sqrt(1 - z[i-1, j] ^ 2)
        }
      }
    }
  )
  
  ### 1
  w1 <- matrix(0, d, d)
  
  system.time(
    
    for(i in 1:d) {
      for(j in 1:d) {
        if((i == 1) & (i == j)) w1[i, j] <- 1
        if((i == 1) & (i < j)) w1[i, j] <- z[i, j]
        if((i > 1) & (i <= j)) w1[i, j] <- z[i, j] * w1[i-1, j] * sqrt(1 - z[i-1, j] ^ 2)
      }
    }
    
  )
  
  
  ### 2
  w2 <- matrix(0, d, d)
  
  system.time(
    
    for(i in 1:d) {
      for(j in 1:d) {
        
        if(i > j) w2[i, j] <- 0
        if((i == j) & (i == 1)) w2[i, j] <- 1
        if((i > 1) & (i == j)) w2[i, j] <- prod(sqrt(1 - z[1:(i-1), j] ^ 2))
        if((i == 1) & (i < j)) w2[i, j] <- z[i, j]
        if((i > 1) & (i < j)) w2[i, j] <- z[i, j] * prod(sqrt(1 - z[1:(i-1), j] ^ 2))
        
      }
    }
    
  )
  
  #### 3
  w3 <- matrix(0, d, d)
  
  system.time(
    {
      w3[1, ] <- c(1, z[1, 2:d])
      
      for(i in 2:d) {
        # for(j in 1:d) {
        for(j in i:d) {
          
          if(i == j) w3[i, j] <- prod(sqrt(1 - z[1:(i-1), j] ^ 2))
          if(i < j) w3[i, j] <- z[i, j] * prod(sqrt(1 - z[1:(i-1), j] ^ 2))
          
        }
      }
    }
  )
  
  
  ### 4
  w4 <- matrix(0, d, d)
  system.time(
    {
      w4[1, ] <- c(1, z[1, 2:d])
      
      for(i in 2:d) {
        for(j in i:d) {
          if(i < j) w4[i, j] <- z[i, j] * w4[i-1, j] * sqrt(1 - z[i-1, j] ^ 2)
          if(i == j) w4[i, j] <- w4[i-1, j] * sqrt(1 - z[i-1, j] ^ 2)
        }
      }
    }
  )
  
  ### compare
  
  w0; w1; w2; w3; w4
  # times
  # 0.007; 0.015; 0.024; 0.015
  (rho0 <- t(w0) %*% w0)
  (rho1 <- t(w1) %*% w1)
  (rho2 <- t(w2) %*% w2)
  (rho3 <- t(w3) %*% w3)
  
  all.equal(rho2, rho3); all.equal(rho0, rho1)
  
  is_positive_definite(rho0); isSymmetric(rho0)
  is_positive_definite(rho1); isSymmetric(rho1)
  is_positive_definite(rho2); isSymmetric(rho2)
  is_positive_definite(rho3); isSymmetric(rho3)
}
