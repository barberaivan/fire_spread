# Evaluate whether it is reasonable to use local-approximate GPs to emulate
# the abc acceptance probability. That would allow to use it as a likelihood
# and easily fit a mixed model. 

library(laGP)
library(randtoolbox)


# Functions ---------------------------------------------------------------

wave_plot <- function(data, response = "abc_prob", alpha = 0.3, best = NULL,
                      x = "par_values", thres = NULL, bin = FALSE,
                      tit = NULL, rc = c(2, 3)) {
  
  if(!is.null(best)) {
    data <- data[order(data[, response], decreasing = TRUE), ]
  }
  
  if(bin) {
    data$bin <- jitter(as.numeric(data$abc_prob >= thres), factor = 0.2)
    response <- "bin"
  }
  
  yy <- range(data[, response])
  yy[2] <- yy[2] * 1.05
  
  # title for first plot
  if(is.null(tit)) {
    tit <- deparse(substitute(data))
  }
  
  par(mfrow = c(rc[1], rc[2]))
  for(i in 1:n_coef) {
    mm <- ifelse(i == 1, tit, NA)
    plot(data[, response] ~ data[, x][, i], ylab = response,
         xlab = par_names[i], ylim = yy, main = mm,
         pch = 19, col = rgb(0, 0, 0, alpha))
    
    if(!is.null(best)) {
      points(data[1:best, response] ~ data[, x][1:best, i],
             pch = 19, col = rgb(0, 0, 1, alpha))
    }
    
    if(!is.null(thres) & response == "abc_prob") {
      abline(h = thres, col = "red")
    }
  }
  par(mfrow = c(1, 1))
}



# tt ----------------------------------------------------------------------


N <- 1e4
ytrue <- 10
yobs <- ytrue + rnorm(N, sd = 0.1)
x <- matrix(seq(0, 1, length.out = N), ncol = 1)

train <- sample(1:N, size = floor(N * 0.9), replace = F)

gp1 <- laGP(x[-train, , drop = F], start = 6, end = 50, 
            X = x[train, , drop = F], Z = yobs[train], method = "mspe")

xout <- matrix(seq(-10, 10+1, length.out = 100))
gp2 <- laGP(xout, start = 6, end = 50, 
            X = x[train, , drop = F], Z = yobs[train], method = "alc")

plot(yobs[train] ~ x[train, ], 
     ylim = c(0, 12),#c(ytrue - 2, ytrue + 2),
     xlim = c(-10, 11),
     col = rgb(0, 0, 0, 0.1), pch = 19)
points(yobs[-train] ~ x[-train, ],
     col = rgb(1, 0, 0, 0.1), pch = 19)
points(gp1$mean ~ x[-train, ],
       col = rgb(0, 0, 1, 0.3), pch = 19)
points(gp2$mean ~ xout[, 1],
       col = rgb(0, 1, 0, 1), pch = 19)

gp2$s2
## Interesting, cuando no hay datos, se va a cero. Entonces, es más seguro
## usarlo en la escala de la prob, porque donde hay pocos datos es porque 
## la proba es generalmente baja. 

# Ahora probamos una función puntuda ---------------------------------

N <- 1e5
xl <- c(-10, 10)
noise_factor <- 0.000001
ss <- diff(xl) * 0.1

x <- matrix(seq(xl[1], xl[2], length.out = N))
mu <- dnorm(x, sd = ss, log = F)
noise <- diff(range(mu)) * noise_factor
y <- rnorm(N, mu, noise)

train <- sample(1:N, size = N * 0.9)

gp1 <- laGP(x[-train, , drop = F], start = 6, end = 50, 
            X = x[train, , drop = F], Z = y[train], method = "alc")

plot(mu ~ x, type = "l", ylim = range(y))
points(y ~ x, col = rgb(0, 0, 0, 0.1), pch = 19, cex = 0.5)
points(gp1$mean ~ x[-train, ],
       col = rgb(0, 0, 1, 0.3), pch = 19, cex = 0.5)


points(y[-train] ~ x[-train, ],
       col = rgb(0, 1, 0, 0.3), pch = 19, cex = 0.5)

points((gp1$mean + 2 * sqrt(gp1$s2)) ~ x[-train, ],
       col = rgb(0, 0, 1, 0.3), pch = 19, cex = 0.5)


## no anda, quizás debería probar con el set de datos de ejemplo, a ver
## si estoy omitiendo algo en el código.


# El problema es que está hecho para simuladores determinísticos, no estocásticos
# Y aún así, si la superficie es muy puntuda, necesita muchos puntos


# Ejemplos paper ----------------------------------------------------------
# Gramacy 2016

x <- seq(-2, 2, by = 0.02)
X <- as.matrix(expand.grid(x, x))
N <- nrow(X)

f2d <- function(x) {
  g <- function(z) return(exp(-(z - 1)^2) + exp(-0.8 * (z + 1)^2) -
                              0.05 * sin(8 * (z + 0.1)))
  -g(x[, 1]) * g(x[, 2])
}
Y <- f2d(X)

Xref <- matrix(c(-1.725, 1.725), nrow = 1)
p.mspe <- laGP(Xref, 6, 50, X, Y, d = 0.1, method = "mspe")
p.alc <- laGP(Xref, 6, 50, X, Y, d = 0.1, method = "alc")

Xi <- rbind(X[p.mspe$Xi, ], X[p.alc$Xi, ])
plot(X[p.mspe$Xi, ], xlab = "x1", ylab = "x2", type = "n",
     main = "comparing local designs", xlim = range(Xi[ ,1]),
     ylim = range(Xi[ ,2]))
text(X[p.mspe$Xi, ], labels = 1:length(p.mspe$Xi), cex = 0.7)
text(X[p.alc$Xi, ], labels = 1:length(p.alc$Xi), cex = 0.7, col = 2)
points(Xref[1], Xref[2], pch = 19, col = 3)
legend("topright", c("mspe", "alc"), text.col = c(1, 2), bty = "n")

p <- rbind(c(p.mspe$mean, p.mspe$s2, p.mspe$df),
           c(p.alc$mean, p.alc$s2, p.alc$df))
colnames(p) <- c("mean", "s2", "df")
rownames(p) <- c("mspe", "alc")
p




xx <- seq(-1.97, 1.95, by = 0.04)
XX <- as.matrix(expand.grid(xx, xx))
YY <- f2d(XX)


nth <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
if(is.na(nth)) nth <- 2
print(nth)

P.alc <- aGP(X, Y, XX, omp.threads = nth, verb = 0)


med <- 0.51
zs <- XX[, 2] == med
sv <- sqrt(P.alc$var[zs])
r <- range(c(-P.alc$mean[zs] + 2 * sv, -P.alc$mean[zs] - 2 * sv))
plot(XX[zs,1], -P.alc$mean[zs], type = "l", lwd = 2, ylim = r,
     xlab = "x1", ylab = "predicted & true response", bty = "n",
     main = "slice through surface")
lines(XX[zs, 1], -P.alc$mean[zs] + 2 * sv, col = 2, lty = 2, lwd = 2)
lines(XX[zs, 1], -P.alc$mean[zs] - 2 * sv, col = 2, lty = 2, lwd = 2)
lines(XX[zs, 1], YY[zs], col = 3, lwd = 2, lty = 3)


# Prueba con el fuego más peque -------------------------------------------

n_coef <- 6
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")
fff <- readRDS("files/overlaps_FI/CoCampana-simulations_list.rds")


dat <- rbind(fff$like_sim,
             fff$abs_prob_check)
dat$ll <- dat$abc_prob ^ 10

wave_plot(dat, response = "ll", alpha = 0.05)

# try laGP
test <- sample(1:nrow(dat), size = nrow(dat) * 0.1)
gg <- laGP(dat$par_values[test, ], 6, 50, 
           dat$par_values[-test, ], dat$ll[-test],
           method = "nn")

plot(gg$mean ~ dat$ll[test])
abline(0, 1, col = 2)
# terrible


# Try NN regression (simplest) --------------------------------------------

# Un KNN? -----------------------------------------------------------------

library(caret)
library(FNN)

knn1 <- knnreg(dat$par_values, dat$ll, 10)
pred_y <- predict(knn1, data.frame(dat$par_values))
plot(pred_y ~ dat$ll); abline(0, 1, col = "red")


# Get KNN
k1 <- get.knn(dat$par_values, k = 10, algorithm = c("kd_tree"))


knn_localreg <- function(xpred, k = 30, weighted = T, repeating = T,
                         X = dat$par_values, Y = dat$ll) {
  
  # ## TESTO
  # xpred <- dat$par_values[25, , drop = F]
  # k = 30
  # X = dat$par_values#[-25, ]
  # Y = dat$ll#[-25]
  # weighted = T
  # repeating = T
  # ##
  
  nn <- get.knnx(X, query = xpred, k = k, algorithm = c("kd_tree"))
  
  if(repeating) {
    nn$nn.index <- nn$nn.index[-1]
    nn$nn.dist <- nn$nn.dist[-1]
  }
  
  xtrain <- X[nn$nn.index, ]
  xtrain2 <- xtrain ^ 2
  xtrain_ <- cbind(xtrain, xtrain2) |> as.data.frame()
  
  xpred_ <- cbind(xpred, xpred ^ 2) |> as.data.frame()
  colnames(xtrain_) <- colnames(xpred_) <- 
    c(par_names, paste(par_names, 2, sep = "_"))
  
  ytrain <- Y[nn$nn.index]
  if(weighted) w = 1 / nn$nn.dist else 1
  
  mm <- lm(ytrain ~ ., data = xtrain_, weights = w)
  out <- predict(mm, newdata = xpred_)
  return(out)
}


yhat <- apply(dat$par_values, 1, function(x) {
  knn_localreg(matrix(x, nrow = 1))
})
# tarda una banda, isnt worth it
