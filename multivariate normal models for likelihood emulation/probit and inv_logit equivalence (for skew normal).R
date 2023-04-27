# approximate Phi with inv_logit
library(mgcv)

N <- 1e4
x <- seq(-30, 30, length.out = N)
mu_real <- exp(pnorm(x, log.p = T))
range(mu_real)

mod <- gam(mu_real ~ x, family = betar(link = "logit"))
coef(mod)
plot(fitted(mod) ~ mu_real); abline(0, 1)

curve(pnorm(x), from = -50, to = 50, n = 1000)
curve(plogis(2 * x), from = -50, to = 50, n = 1000, add = TRUE, col = 2)

xlims <- seq(-500, -20, by = 1)
mm <- cbind(xlims, pnorm(xlims), plogis(xlims * 2))
colnames(mm) <- c("x", "phi", "invlogit")
View(mm)

curve(pnorm(x), from = -5, to = 5, n = 1000)
curve(plogis(1.9 * x), n = 1000, add = TRUE, col = 2)

ff <- function(b) {
  ext <- 10
  xlims <- seq(-ext, ext, length.out = 10000)
  a1 <- plogis(b * xlims)
  a2 <- pnorm(xlims)
  return(sum((a1 - a2) ^ 2))
}

op <- optim(1.9, ff, lower = 1, upper = 10, method = "Brent")
op$par

curve(pnorm(x), from = -5, to = 5, n = 1000)
curve(plogis(op$par * x), n = 1000, add = TRUE, col = 2)

# Discusión interesante acá:
# https://discourse.mc-stan.org/t/normal-lcdf-underflows-much-faster-than-pnorm-log-true/4581/3

# una aprox a log(phi) es
log_phi <- function(x) {
  cc <- log(1 / sqrt(2 * pi)) - 0.5 * x ^ 3 - log(-x)
  return(cc)
}

xlims <- seq(-500, -20, by = 1)
plot(log_phi(xlims) ~ pnorm(xlims, log.p = T))
