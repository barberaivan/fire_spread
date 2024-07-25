
# Refit better GAMs -------------------------------------------------------

# testear 2015_42

dterms <- readxl::read_excel(file.path("data", "focal fires data",
                                       "gam_terms_to_include.xls"),
                             sheet = 2)
dterms[is.na(dterms)] <- 0
# replace NA with zero
all(colnames(dterms)[-1] == names(gam_terms))

simfiles <- list.files(target_dir, pattern = "-simulations_list.rds")



##### TO LOOP BY FIRES

# 2015_42
# 2015_16
# 2015_11
# 2014_1
# 2013_32
wlist <- readRDS(file.path(target_dir, "2013_32-simulations_list.rds"))

# Get previous support

support <- rbind(params_lower, params_upper)
support[2, "steps"] <- max(wlist$like_sim$par_values[, "steps"])



# How many k are needed? --------------------------------------------------

like_sim <- wlist$like_sim
thres <- wlist$thres

above <- like_sim$overlap >= thres

yvar <- "shrubland"
xvar <- "slope"

plot(like_sim$par_values[above, yvar] ~ like_sim$par_values[above, xvar],
     xlim = support[, xvar], ylim = support[, yvar])

convhull <- chull(like_sim$par_values[above, xvar],
                  like_sim$par_values[above, yvar])
temp <- like_sim$par_values[above, c(xvar, yvar)]
borders <- temp[convhull, ]

gg <- expand.grid(
  xvar = seq(support[1, xvar], support[2, xvar], length.out = 150),
  yvar = seq(support[1, yvar], support[2, yvar], length.out = 150)
)

inside_grid <- sp::point.in.polygon(gg$xvar, gg$yvar,
                                    borders[, xvar], borders[, yvar])
gg$y <- inside_grid

# Fit bivariate GAM
kmarg <- 15
kint <- 10

mtest <- bam(
  y ~ s(yvar, k = kmarg, bs = "cr") + s(xvar, k = kmarg, bs = "cr") +
    ti(yvar, xvar, k = kint, bs = "cr"),
  data = gg, method = "fREML", discrete = T, nthreads = 8,
  family = "binomial"
)

gg$p <- fitted(mtest)

ggplot(gg) +
  geom_tile(aes(xvar, yvar, fill = p)) +
  scale_fill_viridis() +
  geom_point(aes(xvar, yvar), data = gg[gg$y == 1, ],
             color = "black", size = 0.5,
             alpha = 0.5)

# Con k = 10 o k = 5 se alcancan cosas muy buenas cuando la forma
# no es muy compleja. No hace falta achicar el espacio.

# Incluso con 2013_32, que tiene la posterior más finita del mundo,
# no hace falta mucho k. Usamos kint = 10 por cuidado, pero con
# kint = 5 hace casi lo mismo.


# lo mismo pero con los datos

gg <- data.frame(
  xvar = like_sim$par_values[, xvar],
  yvar = like_sim$par_values[, yvar],
  y = as.numeric(like_sim$overlap >= thres)
)

ggrid <- expand.grid(
  xvar = seq(support[1, xvar], support[2, xvar], length.out = 150),
  yvar = seq(support[1, yvar], support[2, yvar], length.out = 150)
)

mtest <- bam(
  y ~ s(yvar, k = kmarg, bs = "cr") + s(xvar, k = kmarg, bs = "cr") +
    ti(yvar, xvar, k = kint, bs = "cr"),
  data = gg, method = "fREML", discrete = T, nthreads = 8,
  family = "binomial"
)

ggrid$p <- predict(mtest, ggrid, type = "response")

ggplot(ggrid) +
  geom_tile(aes(xvar, yvar, fill = p)) +
  scale_fill_viridis() +
  geom_point(aes(xvar, yvar), data = gg[gg$y == 1, ],
             color = "black", size = 1,
             alpha = 0.5)





# reduce support and filter data
support_reduced <- reduce_support(
  wlist$like_sim$par_values[wlist$like_sim$overlap >= wlist$thres, ],
  support, prop = 0.5
)

apply(wlist$like_sim$par_values[wlist$like_sim$overlap >= wlist$thres, ],
      2 , range)

rows_inside <- within_support(wlist$like_sim$par_values, support_reduced)
like_sim_red <- wlist$like_sim[rows_inside, ]

above <- like_sim_red$overlap >= wlist$thres
below <- like_sim_red$overlap < wlist$thres
all_below <- wlist$like_sim$overlap < wlist$thres

plot(wlist$like_sim$par_values[, "dry"] ~ wlist$like_sim$par_values[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"])

plot(like_sim_red$par_values[above, "dry"] ~ like_sim_red$par_values[above, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"])

plot(like_sim_red$par_values[, "dry"] ~ like_sim_red$par_values[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"])


plot(wlist$like_sim$par_values[all_below, "dry"] ~
       wlist$like_sim$par_values[all_below, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"])


# fit GAM

tt <- as.matrix(dterms[dterms$fire_id == "2015_42", terms_names])
terms_use <- colnames(tt)[tt == 1]

### OJO, reducir el soporte puede traer problemas porque el espacio no está
### tan bien explorado en términos de interacciones.
### O sea, quizás simplemente hagan falta más k

data_gam <- cbind(y = as.numeric(wlist$like_sim$overlap >= wlist$thres),
                  as.data.frame(wlist$like_sim$par_values))

k_int = 20; k_side = 20
gam_big <- bam(
  make_gam_formula(terms_use), data = data_gam, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)
ff1 <- fitted(gam_big)
ff2 <- fitted(wlist$gam_bern)

par(mfrow = c(2, 1))
plot(ff1 ~ wlist$like_sim$overlap, pch = 19, col = rgb(0, 0, 0, 0.1))
abline(v = wlist$thres, col = "red")

plot(ff2 ~ wlist$like_sim$overlap, pch = 19, col = rgb(0, 0, 0, 0.1))
abline(v = wlist$thres, col = "red")
par(mfrow = c(1, 1))

#still bad with high k.

# Try by evaluating more particles in the reduced support setting.
# reduce support and filter data
support_reduced <- reduce_support(
  wlist$like_sim$par_values[wlist$like_sim$overlap >= wlist$thres, ],
  support, prop = 0.5
)


# 1000 particles in the new support


ss <- sobol(n = 5000, dim = n_coef)
particles <- scale_params(ss, support_reduced)

fire_file <- size_data$file[45]
fire_name <- size_data$fire_id[45]
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids")]
registerDoMC(15)
wlocal <- similarity_simulate_parallel(particles, spread_data)
wlocal$wave <- max(wlist$like_sim$wave) + 1
wlocal$phase <- "sobol2"
colnames(wlocal$par_values) <- par_names

# merge new and ols dimulations in the reduced support
likered2 <- rbind(like_sim_red,
                  wlocal[, names(like_sim_red)])
aa <- likered2$overlap >= wlist$thres
bb <- likered2$overlap < wlist$thres

plot(likered2$par_values[bb, "dry"] ~ likered2$par_values[bb, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))

points(likered2$par_values[aa, "dry"] ~ likered2$par_values[aa, "wind"],
       pch = 19, col = rgb(1, 0, 0, 0.1))

# GAM on reduced support
tt <- as.matrix(dterms[dterms$fire_id == "2015_42", terms_names])
terms_use <- colnames(tt)[tt == 1]
data_gam <- cbind(y = as.numeric(likered2$overlap >= wlist$thres),
                  as.data.frame(likered2$par_values))
k_int = 15; k_side = 15

gam_big <- bam(
  make_gam_formula(terms_use), data = data_gam, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)
ff1 <- fitted(gam_big)
plot(ff1 ~ likered2$overlap, pch = 19, col = rgb(0, 0, 0, 0.1))
abline(v = wlist$thres, col = "red")

length(coef(gam_big))

# queda chequear nuevo GAM
r_gam <- rejection_sample_parallel(200, gam_big, support_reduced,
                                   centre = 0.5,
                                   cores = 10)
draws <- do.call("rbind", r_gam) %>% as.data.frame

# Check GAM
ids <- sample(1:nrow(draws), size = 1000, replace = F)
ppmat <- as.matrix(draws[ids, ])

overlap_check <- similarity_simulate_parallel(particles = ppmat,
                                              fire_data = spread_data)

hist(overlap_check$overlap, xlim = c(0, 1), main = fire_name,
     xlab = "Overlap", ylab = "Relative frequency", freq = F)
abline(v = wlist$thres, col = 2, lwd = 2)
abline(v = mean(overlap_check$overlap), col = 4, lwd = 2)

sum(overlap_check$overlap >= wlist$thres)



# FIT distribution to good particles.

library(kdevine)

samples_prev <- likered2[likered2$overlap >= wlist$thres &
                           likered2$phase %in% c("sobol", "sobol2"),
                         "par_values"]


tested_ok <- overlap_check[overlap_check$overlap >= wlist$thres, ]
distri <- kdevine(tested_ok$par_values,
                  xmin = support_reduced[1, ],
                  xmax = support_reduced[2, ])
samples <- rkdevine(1000, distri)
ddss <- dkdevine(samples, distri)

overlap_check2 <- similarity_simulate_parallel(particles = samples,
                                               fire_data = spread_data)

hist(overlap_check2$overlap, xlim = c(0, 1), main = fire_name,
     xlab = "Overlap", ylab = "Relative frequency", freq = F)
abline(v = wlist$thres, col = 2, lwd = 2)
abline(v = mean(overlap_check2$overlap), col = 4, lwd = 2)

for(i in 1:n_coef) {
  plot(density(samples[, i], adjust = 1.5,
               from = support_reduced[1, i],
               to = support_reduced[2, i]),
       main = par_names[i])
}

colnames(samples) <- par_names

plot(samples[, "dry"] ~ samples[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))


plot(log(ddss) ~ overlap_check2$overlap)


# Esto anda bastante bien!

# Sample GAM with flat probability ----------------------------------------
# to reduce the GAM-introduced error

## ???


# SMC with kdevine? -------------------------------------------------------

# fit density to the particles above threshold.
# reproduce particles wit weight = 1 / dens

par_start <- wlist$like_sim$par_values[wlist$like_sim$overlap >= wlist$thres, ]
dist0 <- kdevine(par_start,
                 xmin = apply(par_start, 2, min),
                 xmax = apply(par_start, 2, max))
d0 <- dkdevine(par_start, dist0)
plot(d0)
min(d0)
hist(log(d0))

w0 <- 1 / normalize(d0)
plot(w0)

# los weights invertidos son malos

# Y si usamos como proposal muestras planas pero sólo donde la
# density de dist0 sea > q10?

low_ll <- quantile(log(d0), prob = 0.01)
min(log(d0))

support_sample <- reduce_support(par_start, support, prop = 0)
dr1 <- rejection_sample_region(100, dist0, support = support_sample, ll_thres = low_ll)


# Tarda mucho si el kde se estima con tantas. Y si usamos las 1000 mejores?

par_start <- par_start[sample(1:nrow(par_start), 1000, replace = F), ]
dist0 <- kdevine(par_start,
                 xmin = apply(par_start, 2, min),
                 xmax = apply(par_start, 2, max))
d0 <- dkdevine(par_start, dist0)
low_ll <- quantile(log(d0), prob = 0.01)
min(log(d0))

support_sample <- reduce_support(par_start, support, prop = 0)
dr1 <- rejection_sample_region(100, dist0, support = support_sample, ll_thres = low_ll)

plot(par_start[, "dry"] ~ par_start[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))

plot(dr1[, "dry"] ~ dr1[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))

# Es muchísimo más rápido si estimamos el KDE con menos muestras.
# Se ve que es una función basada en datos.



par_start <- par_start[1:1000, ]
dist0 <- kdevine(par_start,
                 xmin = apply(par_start, 2, min),
                 xmax = apply(par_start, 2, max))
d0 <- dkdevine(par_start, dist0)
low_ll <- quantile(log(d0), prob = 0.01)
min(log(d0))

support_sample <- reduce_support(par_start, support, prop = 0)
dr1 <- rejection_sample_region(100, dist0, support = support_sample, ll_thres = low_ll)

plot(par_start[, "dry"] ~ par_start[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))

plot(dr1[, "dry"] ~ dr1[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(1, 0, 1, 0.1))



# IDEAS

# ajustar mvn a las partículas above threshold.
# samplear de una uniforme y luego filtrar las que superan el la mínima
# density (o algún valor bajo) en la mvn
# ver mvnfast::

# Luego evaluar las muestras que quedan con el simulador y ajustar un
# kdevine.
# Si pudiera aproximar todo con MVNs, sería hermoso, podría muestrear la
# joint posterior con Stan.

library(sn)

par_start <- wlist$like_sim$par_values[wlist$like_sim$overlap >= wlist$thres, ]

par_start_logit <- as.data.frame(logit_scaled(par_start, support))

mm <- selm(
  cbind(wet, subalpine, dry, shrubland, grassland, slope, wind, steps) ~ 1,
  data = par_start_logit, family = "ST"
) # error con SN
mm

mus <- apply(par_start_logit, 2, mean)
V <- cov(par_start_logit)

xxlogit <- rmvn(5000, mus, V)
samples <- invlogit_scaled(xxlogit, support)

plot(par_start[, "dry"] ~ par_start[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(0, 0, 1, 0.1))

plot(samples[, "dry"] ~ samples[, "wind"],
     xlim = support_reduced[, "wind"], ylim = support_reduced[, "dry"],
     pch = 19, col = rgb(1, 0, 1, 0.1))




# Usando KDE --------------------------------------------------------------

# Idea: con los puntos muestreados, ajustar un KDE plano.
# O sea, un KDE pero en donde yo diga que la density es constante arriba de un umbral.

# Probemos con 1999_28
wlist <- readRDS(file.path(target_dir, "1999_28-simulations_list.rds"))

# Get previous support
like_sim <- wlist$like_sim
thres <- wlist$thres
above <- like_sim$overlap >= thres
support <- rbind(params_lower, params_upper)
support[2, "steps"] <- max(wlist$like_sim$par_values[, "steps"])

yvar <- "grassland" # shrubland
xvar <- "wind"      # slope

plot(like_sim$par_values[above, yvar] ~ like_sim$par_values[above, xvar],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar)

# Fit bivariate density

# use the 1000 best samples
lss <- like_sim[order(like_sim$overlap, decreasing = T), ]
# ss <- sample(1:10000, 1000, replace = F)
aa <- lss$overlap >= thres
lss <- lss[aa, ]

plot(lss$par_values[, yvar] ~ lss$par_values[, xvar],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar)

X <- lss$par_values[, c(xvar, yvar)]
colnames(X) <- c("x", "y")

dfit <- vine(
  X,
  margins_controls = list(
    xmin = support[1, c(xvar, yvar)], # apply(X, 2, min),#
    xmax = support[2, c(xvar, yvar)]  # apply(X, 2, max) #
  ), # Esto ayuda, pero sigue bastante mal.
  copula_controls = list(
    family_set = "nonparametric",
    nonpar_method = "quadratic",
    selcrit = "loglik"
  ),
  cores = 15
)

dfit2 <- kdevine::kdevine(
  X,
  xmin = support[1, c(xvar, yvar)],
  xmax = support[2, c(xvar, yvar)],
  cores = 15
)

# plot samples
X2 <- rvine(10000, dfit)
X3 <- kdevine::rkdevine(10000, dfit2)
colnames(X3) <- colnames(X)

par(mfrow = c(1, 3))

plot(lss$par_values[, yvar] ~ lss$par_values[, xvar],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "data")

plot(X2[, "y"] ~ X2[, "x"],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "vine")

plot(X3[, "y"] ~ X3[, "x"],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "kdevine")

par(mfrow = c(1, 1))


fitted_dens <- dvine(X, dfit, cores = 10)
hist(log(fitted_dens))
dthres <- quantile(log(fitted_dens), 0.0001)

gg <- expand.grid(
  x = seq(support[1, xvar], support[2, xvar], length.out = 500),
  y = seq(support[1, yvar], support[2, yvar], length.out = 500)
)

gg$p <- dvine(as.matrix(gg), dfit)
# gg$p <- kdevine::dkdevine(as.matrix(gg), dfit2)
gg$logp <- log(gg$p)

dobs <- lss$par_values[, c(xvar, yvar)]
colnames(dobs) <- c("x", "y")

ggplot(gg) +
  geom_tile(aes(x, y, fill = logp)) +
  scale_fill_viridis(limits = c(dthres, 0)) +
  geom_point(aes(x, y), data = dobs, alpha = 0.5)

# Si anda mal, no voy a poder usarlo luego.



# usando las mejores 1000 -------------------------------------------------

lss <- like_sim[order(like_sim$overlap, decreasing = T), ]
lss <- lss[1:1000, ]

plot(lss$par_values[, yvar] ~ lss$par_values[, xvar],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar)

X <- lss$par_values[, c(xvar, yvar)]
colnames(X) <- c("x", "y")

dfit <- vine(
  X,
  margins_controls = list(
    xmin = support[1, c(xvar, yvar)], # apply(X, 2, min),#
    xmax = support[2, c(xvar, yvar)]  # apply(X, 2, max) #
  ), # Esto ayuda, pero sigue bastante mal.
  copula_controls = list(
    family_set = "nonparametric",
    nonpar_method = "quadratic",
    selcrit = "loglik"
  ),
  cores = 15
)

dfit2 <- kdevine::kdevine(
  X,
  xmin = support[1, c(xvar, yvar)],
  xmax = support[2, c(xvar, yvar)],
  cores = 15,
  mult_1d = 3
)

# plot samples
X2 <- rvine(10000, dfit)
X3 <- kdevine::rkdevine(10000, dfit2)
colnames(X3) <- colnames(X)

par(mfrow = c(1, 3))

plot(lss$par_values[, yvar] ~ lss$par_values[, xvar],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "data")

plot(X2[, "y"] ~ X2[, "x"],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "vine")

plot(X3[, "y"] ~ X3[, "x"],
     xlim = support[, xvar], ylim = support[, yvar],
     pch = 19, col = rgb(0, 0, 0, 0.2),
     ylab = yvar, xlab = xvar, main = "kdevine")

par(mfrow = c(1, 1))


# fitted_dens <- dvine(X, dfit, cores = 10)
fitted_dens <- dkdevine(X, dfit2)
hist(log(fitted_dens))

gg <- expand.grid(
  x = seq(support[1, xvar], support[2, xvar], length.out = 500),
  y = seq(support[1, yvar], support[2, yvar], length.out = 500)
)

gg <- expand.grid(
  x = seq(5, 12, length.out = 100),
  y = seq(-50, -40, length.out = 100)
)

# gg$p <- dvine(as.matrix(gg), dfit)
gg$p <- kdevine::dkdevine(as.matrix(gg), dfit2)
gg$logp <- log(gg$p)

dobs <- lss$par_values[, c(xvar, yvar)]
colnames(dobs) <- c("x", "y")

dthres <- quantile(log(fitted_dens), 0.0001)

ggplot(gg) +
  geom_tile(aes(x, y, fill = logp)) +
  scale_fill_viridis(limits = c(dthres, max(log(fitted_dens)))) +
  geom_point(aes(x, y), data = dobs, alpha = 0.5)

# Problem: the kdevine fits a very wiggly density, with many gaps that I
# interpret as noise. Perhaps a multivariate skew-normal is better.