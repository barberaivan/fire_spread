# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(rstan)
library(sn)


# Functions ---------------------------------------------------------------

# Function to transform from original (constrained) support to the
# unconstrained one, using scaled logit and log (for steps)
unconstrain <- function(x, support) {
  xun <- x
  
  names_logit <- colnames(x)[colnames(x) != "steps"]
  for(j in names_logit) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }
  
  xun[, "steps"] <- log(x[, "steps"])
  
  return(xun)
}

# The inverse of unconstrain
constrain <- function(xun, support) {
  xc <- xun
  
  names_logit <- colnames(xun)[colnames(xun) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }
  
  xc[, "steps"] <- exp(xun[, "steps"])
  
  return(xc)
}

make_parnames <- function(col, variable = "x", ids = 1:6) {
  # variable = "x"; col = 66; ids = c(1, 2, 3)
  paste(variable, "[", ids, ", ", col, "]", sep = "")
  # 
  # "x[1, 1]", "x[2, 1]", "x[3, 1]", 
  # "x[4, 1]", "x[5, 1]", "x[6, 1]")
}

# Constants ---------------------------------------------------------------

par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")
n_coef <- K <- length(par_names)

# load file with constants to standardize
fi_params <- readRDS(file.path("data", "NDVI_regional_data",
                               "flammability_indices.rds"))
slope_sd <- fi_params$slope_term_sd

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

# Import particles (samples later) ----------------------------------------

data_dir <- file.path("files", "overlaps_FI")
fnames <- list.files(data_dir, pattern = "-simulations_list.rds")
n_fires <- J <- length(fnames)

data_list <- lapply(1:length(fnames), function(f) {
  # f = 1
  dl <- readRDS(file.path(data_dir, fnames[f]))
  return(dl)
})

draws_list <- lapply(1:length(fnames), function(f) {
  dl <- data_list[[f]]
  parts <- rbind(dl$like_sim[, colnames(dl$abc_prob_check)], 
                 dl$abc_prob_check)
  draws <- parts[parts$abc_prob >= dl$thres, "par_values"] 
  return(draws)
})

draws_list_unc <- lapply(1:length(fnames), function(f) {
  support <- rbind(params_lower, params_upper)
  support[2, "steps"] <- max(data_list[[f]]$like_sim$par_values[, "steps"])
  parts <- unconstrain(draws_list[[f]], support)
  
  # remove non-finite values
  ids_ok <- sapply(1:nrow(parts), function(i) {
    return(all(is.finite(parts[i, ])) && !anyNA(parts[i, ]))
  })
  
  return(parts[ids_ok, ])
})

# Estimate MVN densities from samples

mus <- matrix(NA, n_coef, n_fires)
Vs <- array(NA, dim = c(n_coef, n_coef, n_fires))
Ps <- Vs

for(j in 1:n_fires) {
  mus[, j] <- colMeans(draws_list_unc[[j]])
  Vlocal <- cov(draws_list_unc[[j]])
  Vs[, , j] <- Vlocal
  Ps[, , j] <- solve(Vlocal)
  if(anyNA(Vlocal) | anyNA(Ps[, , j]) | anyNA(mus[, j])) print(j)
}

# Model1: fully hierarchical model ----------------------------------------

sdata1 <- list(J = n_fires, K = n_coef, mus = mus, Ps = aperm(Ps, c(3, 1, 2)))

smodel1 <- stan_model("hierarchical model/normal_pseudolikelihoods_01.stan")

m1 <- sampling(
  smodel1, data = sdata1, 
  cores = 4, chains = 4, iter = 2000
)

# Model2: mixed model ----------------------------------------------------

K <- n_coef
ids <- 1:K
idfix <- 1:4
Kf <- length(idfix)
Kr <- K - Kf
idran <- ids[!(ids %in% idfix)]

sdata2 <- list(J = n_fires, K = n_coef, 
               mus = mus, Ps = aperm(Ps, c(3, 1, 2)),
               Kf = Kf, idfix = idfix,
               Kr = Kr, idran = idran)

smodel2 <- stan_model("hierarchical model/normal_pseudolikelihoods_02.stan")

m2 <- sampling(
  smodel2, data = sdata2, 
  cores = 4, chains = 4, iter = 2000
)

# Si sólo hay un random o un fixed, no le gusta el código porque me dice que 
# declaré un array (en idfix o idran) pero le estoy pasando algo de cero dim.



# Fit multivariate skew-normal distributions} -----------------------------

xis <- matrix(NA, n_coef, n_fires)
Omegas <- array(NA, dim = c(n_coef, n_coef, n_fires))
Psn <- Omegas
alphas <- xis
sigmas_inv <- xis

selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

# save dp parameters for each fire in a list, to visualize with sn
sn_list <- vector("list", J)

for(j in 1:n_fires) {
  print(j)
  X <- draws_list_unc[[j]]
  
  # Fit multivariate skew-normal distribution
  snfit <- sn::selm(formula = selm_form,
                    data = as.data.frame(X), family = "SN")
  dpvec <- coef(snfit, "DP", vector = TRUE)
  
  # If it didn't work, use a normal aproximation
  if(anyNA(dpvec) | !all(is.finite(dpvec))) {
    print(paste("Normal approximation in", j))
    
    dp <- list(xi = colMeans(X),
               Omega = cov(X),
               alpha = rep(0, K))
    
    distrib <- makeSECdistr(dp, "SN", "Loglik", par_names)
    sn_list[[j]] <- distrib
    
    xis[, j] <- dp$xi
    Omegas[, , j] <- dp$Omega
    Psn[, , j] <- solve(dp$Omega)
    alphas[, j] <- dp$alpha
    sigmas_inv[, j] <- 1 / sqrt(diag(dp$Omega))
  } else {
    dp <- coef(snfit, "DP", vector = FALSE)
    names(dp)[1] <- "xi"
    dp$xi <- dp$xi[1, ]
    
    distrib <- makeSECdistr(dp, "SN", "Loglik", par_names)
    sn_list[[j]] <- distrib
    
    xis[, j] <- dp$xi
    Omegas[, , j] <- dp$Omega
    Psn[, , j] <- solve(dp$Omega)
    alphas[, j] <- dp$alpha
    sigmas_inv[, j] <- 1 / sqrt(diag(dp$Omega))
  }
}

# Model3: hierarchical skew-normal ----------------------------------------

sdata3 <- list(J = n_fires, K = n_coef, 
               xis = xis, 
               Ps = aperm(Psn, c(3, 1, 2)),
               alphas = alphas,
               sigmas_inv = sigmas_inv)

smodel3 <- stan_model("hierarchical model/normal_pseudolikelihoods_03.stan")

make_inits <- function() {
  xstart <- matrix(NA, K, J)
  for(j in 1:J) {
    xstart[, j] <- rnorm(K, mus[, j], sqrt(diag(Vs[, , j])) * 0.1)
  }
  out <- list(x = xstart)
  return(out)
}

m3 <- sampling(
  smodel3, data = sdata3, thin = 100,
  cores = 1, chains = 1, iter = 15000, warmup = 5000,
  control = list(adapt_delta = 0.95),
  init = make_inits
)

sm3 <- summary(m3)[[1]]
View(sm3)
saveRDS(m3, "files/hierarchical_model/stan_model_m3.rds")
# Warning messages:
#   1: There were 917 divergent transitions after warmup. See
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 33 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems

# Tardó como media hora


pairs(m3, pars = c("x[1, 1]", "x[2, 1]", "x[3, 1]", 
                   "x[4, 1]", "x[5, 1]", "x[6, 1]"))

# El muestreo es muy ineficiente, casi que obtiene una muestra 
# buena cada 4000 / 50 iter 

j <- 30
pairs(m3, pars = make_parnames(j))
plot(sn_list[[j]])
# muchas div transitions, pero muestrea todo muy bien.


# Skew-normal scaling -----------------------------------------------------

pp <- ppoints(100)
sn0 <- qsn(pp, xi = 0, omega = 1, alpha = 3)
sn2 <- qsn(pp, xi = 4, omega = 2, alpha = 3)
sn2z <- (sn2 - 4) / 2

# plot(sn0 ~ sn2z); abline(0, 1, col = 2)

s1 <- rsn(1e4, xi = 0, omega = 1, alpha = 3)
# plot(density(s1))
# plot(density(s1 * 6 - 4))
mean(s1)
sd(s1)

?extractSECdistr()
vcov(sn_list[[1]])
mean(sn_list[[1]])
sn_list[[1]]


# Model4: hierarchical skew-normal scaled ---------------------------------

## Create scaled matrices. 
# for the sn to have mean = 0 and sd = 1, 
# the xi (position) should be subtracted the mean, 
# and the omegas should be divided by the scale.

xis_z <- matrix(NA, n_coef, n_fires)
Omegas_z <- array(NA, dim = c(n_coef, n_coef, n_fires))
Psn_z <- Omegas
alphas <- xis
sigmas_inv_z <- xis

selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

# save dp parameters for each fire in a list, to visualize with sn
sn_list_z <- vector("list", J)

for(j in 1:n_fires) {
  # j = 1
  print(j)
  X <- draws_list_unc[[j]]
  
  # Fit multivariate skew-normal distribution
  snfit <- sn::selm(formula = selm_form,
                    data = as.data.frame(X), family = "SN")
  dpvec <- coef(snfit, "DP", vector = TRUE)
  
  dp <- coef(snfit, "DP", vector = FALSE)
  names(dp)[1] <- "xi"
  dp$xi <- dp$xi[1, ]
  
  distrib <- makeSECdistr(dp, "SN", "Loglik", par_names)
  
  means <- mean(distrib)
  V <- vcov(distrib)
  R <- cov2cor(dp$Omega)
  s <- sqrt(diag(V))
  
  omega_scaled <- sqrt(diag(dp$Omega)) / s
  Omega_scaled <- outer(omega_scaled, omega_scaled) * R
  
  # Force symmetry in Omega_scaled
  Omega_scaled[lower.tri(Omega_scaled)] = t(Omega_scaled)[lower.tri(Omega_scaled)]
  
  xis_z[, j] <- dp$xi - means
  Omegas_z[, , j] <- Omega_scaled
  Psn_z[, , j] <- solve(Omega_scaled)
  alphas[, j] <- dp$alpha
  sigmas_inv_z[, j] <- 1 / omega_scaled
  
  dpz <- list(xi = dp$xi - means, Omega = Omega_scaled, alpha = dp$alpha)
  distrib_z <- makeSECdistr(dpz, "SN", "Loglik_centred", par_names)
  sn_list_z[[j]] <- distrib_z
}


## Stan stuff

sdata4 <- list(J = n_fires, K = n_coef, 
               xis = xis_z, 
               Ps = aperm(Psn_z, c(3, 1, 2)),
               alphas = alphas,
               sigmas_inv = sigmas_inv_z)

# smodel3 <- stan_model("hierarchical model/normal_pseudolikelihoods_03.stan")

make_inits_z <- function() {
  xstart <- matrix(NA, K, J)
  for(j in 1:J) {
    xstart[, j] <- rnorm(K, rep(0, K), 0.1)
  }
  out <- list(x = xstart)
  return(out)
}

m4 <- sampling(
  smodel3, data = sdata4, thin = 100,
  # cores = 4, chains = 4, iter = 2000,
  cores = 1, chains = 1, iter = 50000, warmup = 5000,
  control = list(adapt_delta = 0.95),
  init = make_inits_z
)
# Ahora debería estar muestreando variables con media = 0, sd = 1
# LLeno de transitions, pero no se queja por bajo neff. Y fue más rápido.

j <- 30
pairs(m4, pars = make_parnames(j))
plot(sn_list_z[[j]])

# Parece andar mejor así, pero no sé si podré usar las variables escaladas
# Cuando agregue el modelo mixto.


# Model5: hierarchical skew-normal scaled parameters ----------------------

# Instead of searching a distribution with mean = 0 and sd = 1 on the margins,
# we create distributions with xi = 0 and omega = sigma = sigma_inv = 1.
# I suspect this would not improve sampling per se, but the msn log density
# may be simplified.

# The only parameter we should pass to Stan would be the inverse correlation 
# matrix (Rinv) and the slant vector

Rinv <- array(NA, dim = c(n_coef, n_coef, n_fires))
alphas <- matrix(NA, K, J)

selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

# save dp parameters for each fire in a list, to visualize with sn
sn_list_z <- vector("list", J)

for(j in 1:n_fires) {
  # j = 1
  print(j)
  X <- draws_list_unc[[j]]
  
  # Fit multivariate skew-normal distribution
  snfit <- sn::selm(formula = selm_form,
                    data = as.data.frame(X), family = "SN")
  dpvec <- coef(snfit, "DP", vector = TRUE)
  
  dp <- coef(snfit, "DP", vector = FALSE)
  
  R <- cov2cor(dp$Omega)
  # Force symmetry in R
  R[lower.tri(R)] = t(R)[lower.tri(R)]
  Rinv[, , j] <- solve(R)
  alphas[, j] <- dp$alpha
  
  dpz <- list(xi = rep(0, K), Omega = R, alpha = dp$alpha)
  distrib_z <- makeSECdistr(dpz, "SN", "Loglik_centred", par_names)
  sn_list_z[[j]] <- distrib_z
}
sn_list_z[[1]]@dp$alpha

## Stan stuff

sdata5 <- list(J = n_fires, K = n_coef, 
               Ps = aperm(Rinv, c(3, 1, 2)),
               alphas = alphas)

smodel5 <- stan_model("hierarchical model/normal_pseudolikelihoods_05.stan")

make_inits_z <- function() {
  xstart <- matrix(NA, K, J)
  for(j in 1:J) {
    xstart[, j] <- rmsn(alpha = sn_list_z[[j]]@dp$alpha, 
                        Omega = sn_list_z[[j]]@dp$Omega)
  }
  out <- list(x = xstart)
  return(out)
}

m5 <- sampling(
  smodel5, data = sdata5, thin = 1,
  # cores = 4, chains = 4, iter = 2000,
  cores = 1, chains = 1, 
  iter = 20000, warmup = 5000,
  # iter = 10,
  control = list(adapt_delta = 0.95),
  init = make_inits_z
)
# La otra parametrización era un poquito más rápida, a pesar de que 
# esta hace menos cuentas.

# LLeno de transitions, pero no se queja por bajo neff. 

j <- 52
pairs(m5, pars = make_parnames(j))
plot(sn_list_z[[j]])

# Parece andar mejor así, pero no sé si podré usar las variables escaladas
# Cuando agregue el modelo mixto.


# Model6: mixed skew-normal -----------------------------------------------

# Instead of searching a distribution with mean = 0 and sd = 1 on the margins,
# we create distributions with xi = 0 and omega = sigma = sigma_inv = 1.
# I suspect this would not improve sampling per se, but the msn log density
# may be simplified.

# The only parameter we should pass to Stan would be the inverse correlation 
# matrix (Rinv) and the slant vector

Rinv <- array(NA, dim = c(n_coef, n_coef, n_fires))
alphas <- matrix(NA, K, J)

selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

# save dp parameters for each fire in a list, to visualize with sn
sn_list_z <- vector("list", J)

for(j in 1:n_fires) {
  # j = 1
  print(j)
  X <- draws_list_unc[[j]]
  
  # Fit multivariate skew-normal distribution
  snfit <- sn::selm(formula = selm_form,
                    data = as.data.frame(X), family = "SN")
  dpvec <- coef(snfit, "DP", vector = TRUE)
  
  dp <- coef(snfit, "DP", vector = FALSE)
  
  R <- cov2cor(dp$Omega)
  # Force symmetry in R
  R[lower.tri(R)] = t(R)[lower.tri(R)]
  Rinv[, , j] <- solve(R)
  alphas[, j] <- dp$alpha
  
  dpz <- list(xi = rep(0, K), Omega = R, alpha = dp$alpha)
  distrib_z <- makeSECdistr(dpz, "SN", "Loglik_centred", par_names)
  sn_list_z[[j]] <- distrib_z
}

## Stan stuff

K <- n_coef
ids <- 1:K
idfix <- 2:3
Kf <- length(idfix)
Kr <- K - Kf
idran <- ids[!(ids %in% idfix)]

sdata6 <- list(J = n_fires, K = n_coef, 
               Kf = Kf, idfix = idfix,
               Kr = Kr, idran = idran,
               Ps = aperm(Rinv, c(3, 1, 2)),
               alphas = alphas)

make_inits_z <- function() {
  xstart <- matrix(NA, Kr, J)
  for(j in 1:J) {
    xstart[, j] <- rmsn(alpha = sn_list_z[[j]]@dp$alpha[idran], 
                        Omega = sn_list_z[[j]]@dp$Omega[idran, idran])
  }
  out <- list(x = xstart)
  return(out)
}

smodel6 <- stan_model("hierarchical model/normal_pseudolikelihoods_06.stan")

m6 <- sampling(
  smodel6, data = sdata6, thin = 10,
  cores = 1, chains = 1, 
  iter = 30000, warmup = 10000,
  # iter = 10,
  control = list(adapt_delta = 0.999),
  init = make_inits_z
)
sm6 <- summary(m6)[[1]]
# View(sm6)
# LLeno de transitions, pero no se queja por bajo neff. 
# With high adapt_delta it takes much longer, pero no mejora el muestreo.
# Igual zafa. Con 20000 muestras, el peor neff es ~790, y tampoco tarda tanto.

j <- 3
pairs(m6, pars = make_parnames(j, ids = 1:Kr))
dd <- rmsn(n = 2000, dp = sn_list_z[[j]]@dp)
pairs(dd[, idran], pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8)
pairs(dd[1:10, idran], pch = 19, col = rgb(1, 0, 0, 1), cex = 2, add = T)

# It makes sense if the posterior and the samples differ a bit, because
# in the model, the sampled random effects are conditional on the fixef 
# effects. If some of them were correlated with the now-fixed effects, 
# their distribution might vary. TRUE? Think about a little longer.
# Yes! Consider that now the fixef may have another distribution, not the same
# as in each separate "likelihood".

fixpar <- as.matrix(m6, pars = "xf")
plot(density(fixpar[, 1]))
plot(density(fixpar[, 2]))

# Tareas -------------------------------------------------------------------

# Comparar con más precisión qué forma de parametrizar las msn es más
# rápida/eficiente. Las opciones son como en los modelos 4 y 5:
# 4 = escalamos la distribución para que tenga mean = 0, sd = 1;
# 5 = escalamos los parámetros para achicar las cuentas, pero no implica (4).

# Leer sobre los jakobian adjustments en Stan, para poder escalar todo 
# y también asignar las previas jerárquicas. 

# Siguientes pasos --------------------------------------------------------

# Cuando tenga las muestras del paso 1, evaluar la distribución del overlap
# en ellas y compararla con muestras obtenidas de una aproximación 
# multi-skew-normal.

# Si los overlaps no bajan mucho, y si las formas no son tan distintas,
# go for it.

# Luego ajustar el modelo completamente jerárquico con el Lunn method 
# (mi MCMC) y en Stan, aproximando con msn.
# Comparar las distris posteriores de todo: hiperparams y fixef,
# pero mirando principalmente los fixef, que es lo más relevante.
# Si la cosa no es tan distinta, ir por el modelo mixto usando msn.

# KEMOSION



# density of particles according to msn -----------------------------------

j <- 3
X <- draws_list_unc[[j]]
dens <- dmsn(X, dp = sn_list[[j]]@dp, log = T)
hist(dens, breaks = 30)
abline(v = quantile(dens, c(0.001, 0.005, 0.01)), col = "red")

hist(dens, xlim = c(-42, 6), 
     breaks = seq(-40, 5, by = 5))
summary(dens)
diff(sort(dens))

# probar el percentil 0.5 % en los fuegos chicos, y según cómo vaya la accept prob, 
# seguimos or not.