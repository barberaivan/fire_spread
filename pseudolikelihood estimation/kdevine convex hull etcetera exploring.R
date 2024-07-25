library(mgcv)
library(kdevine)
library(geometry)
library(rvinecopulib) # recomiendan usar este

D <- 8
N <- 1e4
X <- rmvn(n = N, mu = rep(0, D), V = diag(1, D))
dX <- dmvn(t(X), mu = rep(0, D), V = diag(1, D))
hist(exp(dX), breaks = 100)
hist(dX, breaks = 100)
# Es normal, salen muy pocos puntos con alta densidad y muchos con baja.
# Y eso empeora a medida que aumenta la dimensionalidad.


ch <- convhulln(X) # tarda un poco bastante la verdad

X2 <- rmvn(n = N, mu = rep(0, D), V = diag(4, D))
inside <- inhulln(ch, X2) # también tarda mucho, una banda, mucho más que un GAM

str(inside)
sum(inside)

# Si esto no fuera tan lento, podría usarlo para muestrear luego, en vez
# de un GAM


kde <- kdevine(X)
d2 <- dkdevine(X2, kde)
hist(d2, breaks = 100)
plot(kde)
summary(kde)

rkdevine(1000, kde)

kde2 <- vine(X, cores = 2)
rr <- rvine(N, kde2, cores = 2)
dd <- dvine(rr, kde2, cores = 2)
hist(dd)
hist(log(dd))
pairs(rr)
plot(density(rr[, 4]))

# Esto anda hermosoooooooo

# usar todos los puntos arriba del threshold, y muestrear candidatos en base 
# a un umbral de densidad, partiendo de una uniforme.
# El umbral de densidad puede usarse como clasificador de "pasa-nopasa".
# rvinecopulib es mucho más rápido que kdevine.

summary(kde2)

# evaluar el tiempo que toma calcular la density según el nro de muestras usado
# para el ajuste