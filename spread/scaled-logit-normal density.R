library(logitnorm)

logit <- function(x) qlogis(x)
logit_scaled <- function(x, L= 0, U = 1) qlogis((x - L) / (U - L))
unit_scale <- function(x, L = 0, U = 1) (x - L) / (U - L)
unit_scale_d <- function(x, L = 0, U = 1) 1 / (U - L) # for jacobian

logitnorm_pdf <- function(x, mu = 0, sigma = 1, log = F) {
  d <- 
    1 / (sigma * sqrt(2 * pi)) *
    1 / (x * (1 - x)) *
    exp(-(logit(x) - mu) ^ 2 / (2 * sigma ^ 2))  
  # deal with zero densities
  d[is.na(d)] <- 0
  if(log) return(log(d)) else return(d)
}

logitnorm_lpdf <- function(x, mu = 0, sigma = 1) {
  ld <- 
    -log(sigma * sqrt(2 * pi)) +
    -log(x * (1 - x)) +
    -(logit(x) - mu) ^ 2 / (2 * sigma ^ 2)
  # deal with zero densities
  ld[is.na(ld)] <- -Inf
  return(ld)
}

scaled_logitnorm_pdf <- function(x, mu = 0, sigma = 1, L = 0, U = 1, log = F) {
  xs <- unit_scale(x, L, U)
  d <- 
    1 / (sigma * sqrt(2 * pi)) *
    1 / (xs * (1 - xs)) *
    exp(-(logit(xs) - mu) ^ 2 / (2 * sigma ^ 2)) *
    unit_scale_d(x, L, U) ## jacobian adjustment
  # deal with zero densities
  d[is.na(d)] <- 0
  if(log) return(log(d)) else return(d)
}

scaled_logitnorm_lpdf <- function(x, mu = 0, sigma = 1, L = 0, U = 1) {
  xs <- unit_scale(x, L, U)
  ld <- 
    -log(sigma * sqrt(2 * pi)) +
    -log(xs * (1 - xs)) +
    -(logit(xs) - mu) ^ 2 / (2 * sigma ^ 2) +
    log(abs(unit_scale_d(x, L, U))) # log-jacobian adjustment
  # deal with zero densities
  ld[is.na(ld)] <- -Inf
  return(ld)
}

curve(logitnorm_pdf(x, -1), lwd = 2)
curve(dlogitnorm(x, -1), add = T, col = 2, lty = 2)

curve(logitnorm_pdf(x, -1, log = T), lwd = 2)
curve(logitnorm::dlogitnorm(x, -1, log = T), add = T, col = 2, lty = 2)

curve(dlogitnorm(x, -1, log = T), lwd = 2)
curve(logitnorm_lpdf(x, -1), add = T, col = 2, lty = 2)


curve(scaled_logitnorm_lpdf(x, -2, L = 2, U = 5), from = 2, to = 5)
curve(scaled_logitnorm_pdf(x, -2, L = 2, U = 5), from = 2, to = 5)
curve(scaled_logitnorm_pdf(x, -2, L = 0, U = 1), from = 0, to = 1)

integrate(f = function(x) scaled_logitnorm_pdf(x, -2, L = 6, U = 10),
          lower = 6, upper = 10)$value

integrate(f = function(x) scaled_logitnorm_pdf(x, -2, 3, L = 8, U = 10),
          lower = 8, upper = 10)$value