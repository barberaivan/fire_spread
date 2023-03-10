library(tidyverse)
library(randtoolbox)   # sobol sequences
library(DiceKriging)   # Fit Gaussian Processes - 
                       # permite definir el ruido de las simulaciones estocásticas.
library(microbenchmark)

library(foreach)       # parallelization
# library(doMC)

# library(future)      # parallelization again
# library(doFuture)

library(doParallel)
library(doRNG)         # reproducible parallelization

library(Rcpp)
sourceCpp("spread_functions.cpp")


# A few constants ---------------------------------------------------------

d <- 9 # likelihood dimension

# Exploring the sobol sequence ----------------------------------------

m1 <- sobol(n = 1000, dim = d, init = T) # use init = F to grow the seqence.
m2 <- sobol(n = 1000, dim = d, init = F)
m3 <- sobol(n = 1000, dim = d, init = F)

# plot(m1[, 1], m1[, 2])
# points(m2[, 1], m2[, 2], col = "red")
# points(m3[, 1], m3[, 2], col = "blue")

# If the first sample (m1) is initialized, all the following will always be the
# same. i.e., the sequence is fixed. init = F is needed so the following sequences
# are not the same as the first.


# simulate prior function -------------------------------------------------

prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15) { # widened priors
  betas <- c(
    "intercept" = rnorm(1, mu_int, sd_int),   # shrubland logit (reference class)
    "subalpine" = rnorm(1, 0, sd_int),        # veg coefficients
    "wet" = rnorm(1, 0, sd_int),
    "dry" = rnorm(1, 0, sd_int),
    "fwi" = rexp(1, r_z),                     # positive
    "aspect" = rexp(1, r_01),                 # positive (northing)
    "wind" = rexp(1, r_01),                   # positive
    "elevation" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )

  return(betas)
}

prior_sim()

# function to compute prior quantiles from the percentiles, which will be
# obtained as sobols sequences in [0, 1] ^ d
prior_q <- function(p, mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15) {
  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- c("intercept", "subalpine", "wet", "dry",
                   "fwi", "aspect", "wind", "elevation", "slope")

  q[, 1] <- qnorm(p[, 1], mean = mu_int, sd = sd_int)
<<<<<<< Updated upstream

  for(i in 2:4) {
=======
  
  for (i in 2:4) {
>>>>>>> Stashed changes
    q[, i] <- qnorm(p[, i], mean = 0, sd = sd_int)
  }

  names_01 <- which(colnames(q) %in% c("aspect", "wind", "slope"))
  for(i in names_01) { # [0, 1] predictors
    q[, i] <- qexp(p[, i], rate = r_01)
  }

  names_z <- which(colnames(q) %in% c("fwi", "elevation"))
  for(i in names_z) { # standardized predictors
    q[, i] <- qexp(p[, i], rate = r_z)
  }

  return(q)
}


n = 500
p <- sobol(n, d)
w0 <- prior_q(p)

# pairs(w0, pch = 19, col = rgb(0, 0, 0, 0.05))


# GP with DiceKriging -----------------------------------------------------

# # probé varios paquetes antes, pero este resultó ser el más piola para mi caso
# # de simulador estocástico
# 
# n = 20
# p <- sobol(n, d)
# w0 <- prior_q(p)
# 
# colnames(p) <- colnames(w0)
# d0 <- as.data.frame(p)
# d0$y <- rnorm(nrow(d0), 
#               (10) * d0$intercept + (-10) * d0$intercept ^ 2,
#               sd = 1)
# d0$noise <- exp((-3) * d0$intercept + (3) * d0$intercept ^ 2)
# 
# # plot(noise ~ intercept, d0)
# 
# des <- d0[, "intercept", drop = F]
# 
# m2 <- km(y ~ intercept + I(intercept ^ 2), design = des, response = d0$y,
#          nugget.estim = TRUE)
# 
# m4 <- km(y ~ intercept + I(intercept ^ 2), design = des, response = d0$y,
#          noise.var = d0$noise)
# 
# # m3 <- km(~1, design = des, response = d0$y,
# #          nugget.estim = TRUE)
# 
# mod <- m4
# dpred <- data.frame(intercept = seq(-2, 2, length.out = 100))
# preds <- predict(mod, dpred, type = "UK")
# dpred$trend <- as.numeric(preds$trend)
# dpred$mean <- preds$mean
# dpred$lower95 <- preds$lower95
# dpred$upper95 <- preds$upper95
# 
# plot(mean ~ trend, dpred); abline(0, 1) # mean and trend are the same here.
# 
# plot(mean ~ intercept, data = dpred, type = "l", ylim = c(-10, 10))
# points(y ~ intercept, data = d0)
# lines(lower95 ~ intercept, dpred, col = "red")
# lines(upper95 ~ intercept, dpred, col = "red")

# Esto funciona!! permite que las evaluaciones se salgan de rango :)

# deberia usar modelos del tipo m4, y parece que son rapidos.



# Load fires data ---------------------------------------------------------

fire_data <- readRDS(file.path("..", "fire_spread_data", "landscapes_ig-known_non-steppe.rds"))

<<<<<<< Updated upstream
mod <- gam(y ~ s(intercept, k = 6, bs = "gp") +
               s(subalpine, k = 6, bs = "gp") +
               s(wet, k = 6, bs = "gp") +
               s(dry, k = 6, bs = "gp") +
               s(fwi, k = 6, bs = "gp") +
               s(aspect, k = 6, bs = "gp") +
               s(wind, k = 6, bs = "gp") +
               s(elevation, k = 6, bs = "gp") +
               s(slope, k = 6, bs = "gp"),
           data = data, method = "REML")
=======
>>>>>>> Stashed changes

# Computing likelihood for the real data set ------------------------------

<<<<<<< Updated upstream
# It's easier to use a GP directly. We need one that can fit a nugget term, to
# take into account unexplained variability.
# In addition, raw results could be used to account for heteroscedasticity and
# model it. However, that would be more complex.
# Try a first wave and compare the results between 2 gp:
#   including or not a quadratic term for the mean.
# Perhaps the best package by now is GauPro. GPfit is only for deterministic
# simulators. (it would be useful to parameterize FARSITE!)
=======
n_fires <- length(fire_data)
fire_names <- names(fire_data)
>>>>>>> Stashed changes

wind_lay <- which(dimnames(fire_data[[1]]$landscape)$layers == "wind") - 1
elev_lay <- which(dimnames(fire_data[[1]]$landscape)$layers == "elev") - 1

coefs <- c(1000,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)

distances <- rep(30, 8) # sides
distances[c(1, 3, 6, 8)] <- 30 * sqrt(2)

eval_one_particle <- function(particle) {
  
  simil <- array(NA, dim = c(10, 11, n_fires),
                 dimnames = list(
                   replicates = 1:10,
                   metrics = 1:11,
                   fires = fire_names
                 ))
  
  for(f in 1:n_fires) {
    print(paste("simulating fire ", f, ": ", fire_names[f], sep = ""))
    
    simil[, , f] <- emulate_loglik_try(
      landscape = fire_data[[f]]$landscape[, , 1:7], # use the SpatRaster
      burnable = fire_data[[f]]$landscape[, , 8],
      ignition_cells = fire_data[[f]]$ig_rowcol, # 0 indexing
      coef = particle,
      wind_layer = wind_lay,
      elev_layer = elev_lay,
      distances = distances,
      upper_limit = 0.5,
      
      fire_ref = fire_data[[f]],
      n_replicates = 10,
      n_indices = 11
    )
  }
  
  return(simil)
}


# simulate 1 time
part <- c(1000, rep(0, 8)) # terrible coefficients (the worst)

# mbm <- microbenchmark(
#   fun = {ss <- eval_one_particle(part)}, 
#   times = 1
# ) # RAM consumption is really low!
# mbm # ~315.7885 s to evaluate one terrible particle (the most burning one)
# start at ram 3.12 Gb, end at 3.38 Gb

# str(ss)
# range(ss[, 1, ]); hist(ss[, 1, ])
# apply(ss, 2, range) # all reasonable values


# worst scenario:
# 900 * 315.7885 / 3600 = 78.94712 # h for 900 particles, without parallelization.



# Parallelization test ----------------------------------------------------

<<<<<<< Updated upstream
ggplot(table_plot, aes(x = intercept, y = y)) +
  geom_point()

ggplot(table_plot, aes(x = intercept, y = y_hat)) +
  geom_point()
=======
# function to loop across particles
loglik_fire_particle <- function(fire, particle) {
    
    res <- emulate_loglik_try(
      landscape = fire$landscape[, , 1:7], # use the SpatRaster
      burnable = fire$landscape[, , 8],
      ignition_cells = fire$ig_rowcol, # 0 indexing
      coef = particle,
      wind_layer = wind_lay,
      elev_layer = elev_lay,
      distances = distances,
      upper_limit = 0.5,
      
      fire_ref = fire,
      n_replicates = 10,
      n_indices = 11
  )
  
  return(res)
}

part <- c(1000, rep(0, 8)) # terrible coefficients (the worst)
loglik_fire_particle(fire_data[[1]], part)
>>>>>>> Stashed changes

# make particles
n = 20; d = 9
p <- sobol(n, d)
w0 <- prior_q(p)
part_list <- lapply(1:nrow(w0), function(par) w0[par, ])


# Compare 2 functions, sequential and parallel

loglik_eval_seq <- function(fire, particles) {
  # particles must be a matrix! (replicates in rows)
  result <- vector(mode = "list", length = nrow(particles))
  for(part in 1:nrow(particles)) {
    result[[part]] <- loglik_fire_particle(fire, particles[part, ])
  }
  return(result)
}

loglik_eval_par <- function(fire, particles_list) {
  # particles must be a list!
  result <- foreach(pp = particles_list) %dopar% {
    loglik_fire_particle(fire, pp)
  }
  return(result)
}

# benchmark
npar = 10; d = 9
p <- sobol(npar, d)
w0 <- prior_q(p)
w0_list <- lapply(1:nrow(w0), function(par) w0[par, ])

# registerDoMC(16)
# mbm_par <- microbenchmark(
#   sequential = {try_seq <- loglik_eval_seq(fire_data[["2015_50"]], w0)},
#   parallel = {try_par <- loglik_eval_par(fire_data[["2015_50"]], w0_list)},
#   times = 1
# ) 
# mbm_par
# con 100 particulas, par / seq fue 42 / 250
# 80 / 322 para 10 particulas (* 10 replicas) en cholila = 25 %

gc()

# Does RAM crash with 100 particles?
# I dont know if it crashed, but R got blank and the RAM at top.

# npar = 100; d = 9 
# p <- sobol(npar, d)
# w0 <- prior_q(p)
# w0_list <- lapply(1:nrow(w0), function(par) w0[par, ])
# try_par <- loglik_eval_par(fire_data[["2015_50"]], w0_list)

# De todos modos, esto no suena nada a shared memory. Incluso con forking dicen
# que a cada proceso se le manda una copia de los objetos requeridos. 



# monitoring RAM usage ----------------------------------------------------

npar = 100; d = 9 
p <- sobol(npar, d)
w0 <- prior_q(p)
w0_list <- lapply(1:nrow(w0), function(par) w0[par, ])

# compare ram usage for 50 or 200 particles
w50 <- prior_q(sobol(100, d))
w200 <- prior_q(sobol(1000, d))

list_mat <- function(mat) lapply(1:nrow(mat), function(row) mat[row, ])

l50 <- list_mat(w50)
l200 <- list_mat(w200)

# registerDoMC(16); loglik_eval_par(fire_data[[1]], l50)
# registerDoMC(16); loglik_eval_par(fire_data_1, l200)

# peakRAM can't measure parallel processes. Evaluating 2000 particles uses just
# a little bit more RAM than evaluating 1000.


# registerDoMC(16); loglik_eval_par(fire_data[["2015_50"]], l50)

# Esto mismo en R base anduvo re bien, sin crashear, y el peak ram no fue tan
# alto. Pero veamos qué pasa si uso más de 20 partículas.

# Probé con 100 y pasó algo hermoso: anduvo bien, con la RAM a tope todo el tiempo
# pero sin crashear. Pero al final me tiró este warning:

# Warning message:
# In mclapply(argsList, FUN, mc.preschedule = preschedule, mc.set.seed = set.seed,  :
#   scheduled cores 12, 15 did not deliver results, all values of the jobs will be 
# affected

<<<<<<< Updated upstream
d0$y <- rnorm(nrow(w0),
              (10) * d0$intercept + (-10) * d0$intercept ^ 2,
              sd = 1)
=======
# más o menos 15 / 100 resultados fueron NULL. Lo bueno es que uno puede saber en 
# cuáles hubo problemas y volver a evaluar esas partículas.
# esos 15 null, si conté bien, coinciden con "15 did not deliver results".
>>>>>>> Stashed changes

# Por lo que leí, parece que pasa con RAM issues. 

# No estoy benchmarqueando, pero si 85 partículas corrieron en ~15 min (quizás menos),
# una wave de 900 seguro tarda menos que 72 / 4 h. Difire mucho entre partículas 
# muestreadas de la previa y la peor peor partícula que había benchmarqueado antes.


# registerDoMC(12)
# pruebita <- loglik_eval_par(fire_data[["2015_50"]], l50)

# Con 12 cores no llega al tope de RAM, y no tiró errores.
# El timing es hermoso, incluso usando 12 cores. 


# ----------------------------------------------------------------------

# Probar con future (multicore) en R base. Quizás ande mejor que parallel.

# Aunque parece que doFuture es lo mismo que %dopar%, con doMC:
# https://cran.r-project.org/web/packages/doFuture/vignettes/doFuture.html

# Definitivamente hay que aproximar la likelihood fuego por fuego. 
# Agarrar los fuegos más grandes y evaluar con cuál se llega al tope de RAM. 
# Para esos grandes, bajar el n_cores a 12 o 14 (cholila anduvo bien con 12, pero 
# quizás se banca un par más). 

# También puede ser útil no cargar todos los paisajes a la vez sino tenerlos 
# guardados por separado, entonces la función que simula hace esto:

# fire <- readRDS("focal_fire_filename.rds")
# loglik_list <- loglik_eval_par(fire, particle_list)
# saveRDS(loglik_lis, "loglik_w1_firename.rds")
# rm(fire); rm(loglik_list); gc() # libera RAM!

# Aunque sería mejor hacer todas las waves con ese fuego y luego pasar a otro.
# Altas ganas de estimar el modelo derecho y no pelotudear con fuegos simulados.


# Using doRNG for reproducible parallel RNG

cl <- makeCluster(16, type = "FORK")
registerDoParallel(cl)
registerDoRNG(seed = 123)
# foreach(i=1:4, .combine = 'c') %dopar% {rnorm(1, mean = 0, sd = 1)}
# stopImplicitCluster()

# registerDoFuture()
# plan(multicore, workers = 12)
# # the function is the same, but if doFuture is loaded, the future's 
# # %dopar% is used, and not the doMC.
# 
# # allow large export for future (if not, throws an error).
# # https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize
# desired_size_mb <- 1000
# options(future.globals.maxSize = desired_size_mb * 1024 ^ 2)

npar = 100; d = 9 
p <- sobol(npar, d)
w0 <- prior_q(p)
w0_list <- list_mat(w0)

<<<<<<< Updated upstream
colnames(p) <- colnames(w0)
d0 <- as.data.frame(p)
d0$y <- rnorm(nrow(d0),
              (10) * d0$intercept + (-10) * d0$intercept ^ 2,
              sd = 1)
d0$noise <- exp((-3) * d0$intercept + (3) * d0$intercept ^ 2)
=======
bench_16cores <- microbenchmark(
  "a" = {pruebita <- loglik_eval_par(fire_data[["2015_53"]], w0_list)},
  times = 1
); bench_16cores
>>>>>>> Stashed changes


# Timings para simular fuegos en paisajes cualquiera.
# En busca de algo que no tarde tanto. 


# 1999_27j_S tarda 58 s para 100 partículas. una wave de 1000 tardaría 10 min.
# * 15 fuegos simulados... 150 min por wave. 
# Feo... busquemos uno más chico.

# 2015_47N tarda 38 s para 100 partículas
38 * 10 * 10 * 15 / 3600 # = 15.83 h  
# costo de 100 partículas * 10 veces para hacer 1000 * 10 waves * 15 fuegos

# 2015_53 (Alerces 2015)
13.86646 * 10 * 10 * 15 / 3600 # = 5.777692 h  





# prueba para ver qué pasa con doRNG cuando la RAM se va a tope.
# Tira un error y sigue corriendo, igual que future.


# tiró un error porque algún process no pudo dar un result. Y si bien siguió corriendo,
# no guardó ningún resultado. En este sentido es peor que doMC. Ese al menos te 
# guarda un NULL.
# Ahora probamos con 12 cores y con previas ensanchadas.
# Nuevamente, parece que no hace tope de RAM.

# ---------------------------------------------------------------------------
# Warning messages:
# 1: UNRELIABLE VALUE: One of the foreach() iterations (‘doFuture-1’) unexpectedly generated random numbers without declaring so. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, use '%dorng%' from the 'doRNG' package instead of '%dopar%'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, set option 'doFuture.rng.onMisuse' to "ignore".
# -----------------------------------------------------------------------------

# CHEQUEAR ESTO!!!

# pero
sapply(pruebita, function(x) is.null(x)) %>% sum()
# tira 0, o sea que parece que anduvo bien

bench100_12cores # 545.1716 seconds. Es sospechoso!
# 545.1716 * 10 / 3600 = 1.514366 h
# mañana revisar bien ese result

hist(sapply(pruebita, function(x) x[, 1])) # makes sense!


# miramos un poco el overlap en func del intercept y el FWI
sss <- sapply(pruebita, function(x) mean(x[, 1]))

loglik_df <- cbind(
	as.data.frame(w0),
	data.frame(lik = sss, log_lik = log(sss))
)

ggplot(loglik_df, aes(intercept, fwi, color = log_lik)) + 
	geom_point(size = 4) + 
	scale_color_viridis(end = 0.8) + 
	theme_bw()

ggplot(loglik_df, aes(intercept, fwi, color = lik)) + 
	geom_point(size = 4) + 
	scale_color_viridis(end = 0.8) + 
	theme_bw()

plot(lik ~ intercept, loglik_df)     # acá se ve bien el efecto del tamaño del paisaje.
plot(log_lik ~ intercept, loglik_df) # intercepts muy altos son malos, pero no llegan a
# verse tan malos porque el paisaje no les alcanza. Tener cuidado con eso, y chequear 
# si no sesga la estimación. (hay que ver cómo evoluciona la estimación del modelo
# al aumentar las waves). Se me ocurre que la curva de likelihood se podría estimar 
# sin esos puntos que se van al diablo. Con la función cuadrática para la media, 
# eso quizás sería bastante sano. Pero si vamos a estar eliminando los puntos en que 
# el intercept se va al diablo en el fuego más grande, significa que ya sabemos que 
# la likelihood hacia allá arriba solo disminuye. Entonces, mejor no usar partículas
# en la zona en que la likelihood se hace plana en 0.05 (overlap_sp scale).

plot(log_lik ~ intercept, loglik_df)
mm <- lm(log_lik ~ intercept + I(intercept ^ 2), data = loglik_df)
curve(coef(mm)[1] + coef(mm)[2] * x + coef(mm)[3] * x ^ 2, add = TRUE)

# claaaro, tira todo lo bueno para intercepts altísimos!
# Pero bueno, esto es un modelo marginal. Hay que ver qué pasa si se tienen en cuenta 
# todos los params. Una buena puede ser ajustar GP marginales para tantear un poco.
# Igual es muy probable que todo esto se corrija mucho al hacer más waves.
# Como que primero se descartan los intercepts muy bajos, lo cual no está mal, 
# luego deberían descartarse los más altos.

# Mañana puedo probar esto mismo con un fuego más pequeño, para hacer muchas waves
# sin que tarde tanto.

# Sobre RAM:

# El uso de RAM es interesante. Se ve que el sistema la mantiene o la libera según
# necesidad. Luego de esa prueba fallada, se quedó con la RAM re cargada aunque no
# hiciera nada, pero luego la descartó cuando le volví a pedir que laburara mucho. 

# Sobre previas:

# Puse previas muy amplias en la pruebita, porque Pájaro mostró en su paper que el
# beta de forest andaba por ~[-20, -5], con un intercept (matorral) de ~[-10, 0].
# De todos modos, para un fuego único, es probable que ese intercept llegue a valores
# más razonables. Probar un poco con los fuegos más pequeños y con los más grandes
# qué regiones de la likelihood tienen valores más altos. Y en base a eso podemos 
# definir las previas. Aunque la verdad que quizás es más seguro hacer un par de 
# waves extra y listo.

