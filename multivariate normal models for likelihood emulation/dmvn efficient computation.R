library(sn)
library(mgcv)

getAnywhere(dmvn)
dmvn_mgcv <- function (x, mu, V, R = NULL) 
{
  if (is.null(R)) 
    R <- chol(V)
  z <- forwardsolve(t(R), x - mu)
  -sum(log(diag(R))) - log(2 * pi) * length(mu)/2 - if (is.matrix(z)) 
    colSums(z^2)/2
  else sum(z^2)/2
}


getAnywhere(dmsn)
dmsn_sn <- 
function (x, xi = rep(0, length(alpha)), Omega, alpha, tau = 0, 
          dp = NULL, log = FALSE) 
{
  if (!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
    stop("You cannot set both component parameters and dp")
  if (!is.null(dp)) {
    if (length(dp) < 3) 
      stop("wrong length of non-null 'dp'")
    xi <- drop(dp[[1]])
    Omega <- dp[[2]]
    alpha <- dp[[3]]
    tau <- if (length(dp) == 4) 
      dp[[4]]
    else 0
  }
  if (any(abs(alpha) == Inf)) 
    stop("Inf's in alpha are not allowed")
  d <- length(alpha)
  Omega <- matrix(Omega, d, d)
  invOmega <- pd.solve(Omega, silent = TRUE, log.det = TRUE)
  if (is.null(invOmega)) 
    stop("Omega matrix is not positive definite")
  logDet <- attr(invOmega, "log.det")
  x <- if (is.vector(x)) 
    matrix(x, 1, d)
  else data.matrix(x)
  if (is.vector(xi)) 
    xi <- outer(rep(1, nrow(x)), as.vector(matrix(xi, 1, 
                                                  d)))
  if (tau == 0) {
    log.const <- logb(2)
    alpha0 <- 0
  }
  else {
    log.const <- -pnorm(tau, log.p = TRUE)
    O.alpha <- cov2cor(Omega) %*% alpha
    alpha0 <- tau * sqrt(1 + sum(alpha * O.alpha))
  }
  X <- t(x - xi)
  Q <- colSums((invOmega %*% X) * X)
  L <- alpha0 + as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
  logPDF <- (log.const - 0.5 * Q + pnorm(L, log.p = TRUE) - 
               0.5 * (d * logb(2 * pi) + logDet))
  if (log) 
    logPDF
  else exp(logPDF)
}



# Q <- colSums((invOmega %*% X) * X)
# this is the important step.

# Perhaps it's better to parameterize the mvn from the precision because it
# doesn't have to be inverted. 
# However, we would still need the sigma vector for the skew-normal...

# https://search.r-project.org/CRAN/refmans/LaplacesDemon/html/dist.Multivariate.Normal.Precision.Cholesky.html
