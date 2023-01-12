# Compare discrepancy metrics using only simulated fires, based on cholila landscape.

library(terra)
library(tidyverse)
library(Rcpp)
library(viridis)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")
sourceCpp("discrepancy_functions.cpp")
# sourceCpp("discrepancy_functions2.cpp") # mismo error

# Por un bug no funciona la overlap_spvd. Así que la haré a mano.

# import landscape
data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"
land <- readRDS(paste(data_dir, "data_cholila_landscape.R", sep = ""))

# get raster for ncol and nrow
land_raster <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
land_raster

# distances for 30 m resolution
distances <- rep(res(land_raster)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(land_raster)[1] * sqrt(2)

# coefficients 
coefs <- c(100, 
           0, 0, 0, 
           0, 0, 0, 0, 0)

ig_location <- cellFromRowCol(land_raster, 
                              nrow(land_raster) / 2, ncol(land_raster) / 2)

# see fire and burnable (it's important to use the burnable layer to 
# avoid huge fires). Note that cholila was mainly contained by non-burnable
# landscape.
# values(land_raster) <- land[, "burnable"]
# plot(land_raster, col = c("black", "green"))
# values(land_raster)[land[, "burned"] == 1] <- 2
# plot(land_raster, col = c("black", "green", "red")) # burnable

# check ig_location falls in burnable area
land[ig_location, "burnable"] == 1

# cholila size
# (land[, "burned"] %>% sum()) / nrow(land) * 100 
# 3.047275 % del paisaje

# Simulate fires ----------------------------------------------------------

# matrix to be filled with fires
nsim <- 100

# fire_sim <- matrix(NA, nrow = nrow(land), ncol = nsim)
# for(i in 1:100) {
#   print(i)
#   fire_sim[, i] <- simulate_fire_cpp(
#     landscape = land[, 1:7],      # to cpp we pass the values, not the raster
#                                   # and remove burnable and burned layers
#     burnable = land[, "burnable"],
#     ignition_cells = ig_location - 1,
#     n_rowcol = c(nrow(land_raster), ncol(land_raster)),
#     coef = coefs,
#     wind_column = 6 - 1,
#     elev_column = 7 - 1,
#     distances = distances,
#     upper_limit = 0.28 # 0.3 para los 1ros 5, luego 0.20, luego 0.25
#     # el timing cambia muchísimo al pasar de 0.3 a 0.2.
#     # con 0.3 ya se hacen gigantes, pero con 0.2 son re peques
#     # 0.25 también hace fuegos muy muy pequeños
#   )
#   
#   # save by parts every ten iterations
#   if((i %% 10) == 0) saveRDS(fire_sim, "fire_sim_028_temp.R")
# }
# saveRDS(fire_sim, "fire_sim_028.R")
# # done

fire_sim <- readRDS("fire_sim_028.R")


# check sizes
sizes_pix <- colSums(fire_sim)
sizes <- sizes_pix / nrow(land) * 100 # size as landscape area percentaje
hist(sizes, breaks = 30, xlab = "size in percentaje (%)")
# son re pequeños con 0.25!!!
# con 0.28 hay muy peques y muy grandes.

largest_id <- which.max(sizes)
 
# # plot the largest and cholila
# fire_choli <- land_raster
# values(fire_choli) <- land[, "burnable"]
# values(fire_choli)[land[, "burned"] == 1] <- 2
# plot(fire_choli, col = c("black", "green", "red"))
# # --
# fire_big <- land_raster
# values(fire_big) <- land[, "burnable"]
# values(fire_big)[fire_sim[, largest_id] == 1] <- 2
# plot(fire_big, col = c("black", "green", "red")) # quemó casi todo (p = 0.28)


# Compare discrepancies ---------------------------------------------------

# fire_sim.R has very small fires, created with p = 0.25.
# fire_sim_03.R will use p = 0.3, and we hope that makes them larger.

# create array to fill with discrepancies
disc_arr <- array(NA, dim = c(nsim, nsim, 5),
                  dimnames = list(fires = 1:nsim,
                                  fires = 1:nsim,
                                  disc = c("ov_sp", "ov_vd", "d_m_raw", 
                                           "mean_size", "dif_size")))

# for(c in 1:(nsim-1)) {  ## this loop also takes long
#   # print(c)
#   for(r in (c+1):nsim) {
#     print(paste("col:", c, "/ row:", r))
#     disc_arr[r, c, "ov_sp"] <- overlap_sp(fire_sim[, c], fire_sim[, r])
#     disc_arr[r, c, "ov_vd"] <- overlap_vd(fire_sim[, c], fire_sim[, r],
#                                           landscape = land[, 1:7])
#     disc_arr[r, c, "d_m_raw"] <- delta_m_raw(fire_sim[, c], fire_sim[, r],
#                                              landscape = land[, 1:7])
#     disc_arr[r, c, "mean_size"] <- mean(sizes[c], sizes[r])
#     disc_arr[r, c, "dif_size"] <- abs(sizes[c] - sizes[r])
#   }
# }
# saveRDS(disc_arr, "fire_sim_disc_arr_028.R")
disc_arr <- readRDS("fire_sim_disc_arr_028.R")
# str(disc_arr)

clean_arr <- function(x) na.omit(as.numeric(x))
discrep <- apply(disc_arr, 3, clean_arr) %>% as.data.frame()

# compute other discrepancies
discrep$ov_spvd <- discrep$ov_sp * 0.8 + discrep$ov_vd * 0.2
discrep$delta_m <- (1 - discrep$ov_sp) + discrep$d_m_raw

# compare
pairs(discrep, pch = 19, col = rgb(0, 0, 0, 0.05))
GGally::ggpairs(discrep, aes(alpha = 0.1))

# check contribution of overlap_sp into delta_m (probably, very small)
discrep$ov_importance <- (1 - discrep$ov_sp) / discrep$delta_m
hist(discrep$ov_importance, breaks = 10) # mostly low
# is the importance related to the difference? I gues it is
plot(ov_importance ~ delta_m, data = discrep)


# So, as the overlaps are bounded, when combined with delta_burned_area_by_veg,
# they have small importance, mainly when the fires are large. 
ggplot(discrep, aes(x = d_m_raw, y = ov_sp, colour = dif_size)) +
  geom_point(alpha = 0.5, size = 2) + 
  scale_color_viridis(option = "B")
# note that for similar overlap values, d_m_raw may be large or small.

# And overlap, for me, is better than ov_vd, because many fires may have 
# ov_vd = 1 with varying overlap. However, this happens mostly at small fires.
ggplot(discrep, aes(x = ov_vd, y = ov_sp, colour = dif_size)) +
  geom_point(alpha = 0.5, size = 2) + 
  scale_color_viridis(option = "B")


GGally::ggpairs(discrep[, c("ov_sp", "ov_vd", "d_m_raw", "dif_size")], 
                aes(alpha = 0.2))
# ov_sp is clearly more related to size difference than ov_vd.
# d_m_raw is also related to dif_size, but maybe not as much as ov_sp.


hist(discrep$ov_sp)


# TO DO -------------------------------------------------------------------

# repeat everything using the larger fires (p = 0.3)
# (they take soooo long to run!!!)

# Y también repetir lo siguiente

# Average overlap ~ sample size -------------------------------------------

# evaluate how the variance in the mean overlap changes as a function of the
# number of fires used to evaluate it. 
# use seq(10, 100, by = 10) s.

# For each fire and sample size, compute the overlap_mean using 10 bootstrap 
# replicates. Then, compute variances by fire, ss and replicate.

ss <- c(1, 5, seq(10, 100, by = 10))
nsize <- length(ss)
nboot <- 100

ov_avg <- array(NA, dim = c(nsize, nboot, nsim),
                dimnames = list(
                  sample_size = ss,
                  boot_rep = 1:nboot,
                  fire = 1:nsim
                ))

for(f in 1:nsim) { # loop over fires
  # f = 1
  if(sizes_pix[f] > 1) {
    # get all overlaps
    ovs_cols <- disc_arr[, f, "ov_sp"]
    ovs_rows <- disc_arr[f, , "ov_sp"]
    ovs <- na.omit(c(ovs_cols, ovs_rows)) %>% as.numeric
    
    if(length(ovs) != (nsim-1)) stop("check overlap size")
    
    for(s in 1:length(ss)) {
      for(b in 1:nboot) {
        sample_ids <- sample(1:length(ovs), size = ss[s], replace = TRUE)
        ov_avg[s, b, f] <- mean(ovs[sample_ids])
      }
    }
  }
}

# Compute variance across bootstrap replicates

ov_var <- apply(ov_avg, c(1, 3), var, na.rm = TRUE) %>% as.data.frame.table()
ov_var <- ov_var[complete.cases(ov_var), ]
names(ov_var)[3] <- "var"
ov_var$sd <- sqrt(ov_var$var)

ggplot(ov_var, aes(x = sample_size, y = var, color = fire, group = fire)) +
  geom_line()

ov_var_agg <- aggregate(cbind(var, sd) ~ sample_size, ov_var, mean)
ov_var_agg$sample_size <- as.numeric(as.character(ov_var_agg$sample_size))

ggplot(ov_var_agg, aes(x = sample_size, y = var)) +
  geom_line() + 
  geom_point(size = 3)

ggplot(ov_var_agg, aes(x = sample_size, y = sd)) +
  geom_line() + 
  geom_point(size = 3)

# 50 o 40 fuegos serían lo ideal, pero ya entre 5 y 10 cambia mucho


# comments ----------------------------------------------------------------

# Es muy probable que para paralelizar tengamos problemas de RAM. 

# Con p = 0.28:
# Los patrones son similares, aunque es difícil ver porque la distrib de overlap
# es muy disyunta, y la de tamaños también.


# Landscape size and overlap (comparing cholila fire) ---------------------

fire_sim_1 <- readRDS("fire_sim_028.R")
fire_sim_2 <- readRDS("fire_sim.R")
fire_sim <- cbind(fire_sim_1, fire_sim_2)
nsim <- ncol(fire_sim)

sizes_pix <- colSums(fire_sim)
sizes <- sizes_pix / nrow(fire_sim) * 100

# # plot the largest and cholila
# fire_choli <- land_raster
# values(fire_choli) <- land[, "burnable"]
# values(fire_choli)[land[, "burned"] == 1] <- 2
# plot(fire_choli, col = c("black", "green", "red"))
# # --
# fire_big <- land_raster
# values(fire_big) <- land[, "burnable"]
# values(fire_big)[fire_sim[, largest_id] == 1] <- 2
# plot(fire_big, col = c("black", "green", "red")) # quemó casi todo (p = 0.28)

# size cholila is
# 3.047275 % del paisaje

# the largest here (almost all the landscape)
# sizes[largest_id]
# 52.20511

# the % burnable is
# sum(land[, "burnable"]) / nrow(land) * 100
# 68.25776 %

# fire_largest <- land[, "burnable", drop = F]

# overlap cholila and the largest simulated:
# (ov_choli_large <- overlap_sp(land[, "burned"], fire_sim[, largest_id]))

# overlap cholila and largest possible:
# (ov_choli_largest <- overlap_sp(land[, "burned"], land[, "burnable"]))

# compute overlap between cholila and all the simulated fires with p = 0.28.
disc_choli <- matrix(NA, nsim, 5)
colnames(disc_choli) <- c("ov_sp", "ov_vd", "d_m_raw", "size", "dif_size")
disc_choli[, "size"] <- sizes
choli_size <- (land[, "burned"] %>% sum()) / nrow(land) * 100 
disc_choli[, "dif_size"] <- abs(choli_size - sizes)

for(i in 1:nsim) {
  print(i)
  disc_choli[i, "ov_sp"] <- overlap_sp(land[, "burned"], fire_sim[, i])
  disc_choli[i, "ov_vd"] <- overlap_vd(land[, "burned"], fire_sim[, i], land)
  disc_choli[i, "d_m_raw"] <- delta_m_raw(land[, "burned"], fire_sim[, i], land)
}

disc_choli <- as.data.frame(disc_choli)
disc_choli$sizes_pix <- sizes_pix
ov_choli_ecdf <- ecdf(disc_choli$ov_sp)
disc_choli$p <- ov_choli_ecdf(disc_choli$ov_sp)

# get delta_m
disc_choli$delta_m <- (1 - disc_choli$ov_sp + disc_choli$d_m_raw)

# compute size quotient with choli in numerator
disc_choli$size_q_num <- choli_size / sizes

# compute size quotient with choli in denominator
disc_choli$size_q_den <- sizes / choli_size 

# add burn_p value
disc_choli$burn_p <- rep(c(0.28, 0.25), each = 100)

# plots

plot(ov_sp ~ dif_size, disc_choli[disc_choli$sizes_pix, ])
plot(ov_sp ~ size, disc_choli); abline(v = choli_size, col = "red")

plot(delta_m ~ size, disc_choli); abline(v = choli_size, col = "red")

# sizes quotients
par(mfrow = c(1, 2))
plot(size_q_num ~ size, disc_choli, main = "Choli numerator")
abline(v = choli_size, col = "red")
plot(size_q_den ~ size, disc_choli, main = "Choli denominator")
abline(v = choli_size, col = "red")
par(mfrow = c(1, 1))


ggplot(disc_choli, aes(x = size, y = ov_sp, colour = as.factor(burn_p),
                       size = as.factor(burn_p))) +
  geom_point() + 
  geom_vline(xintercept = choli_size)


aggregate(ov_sp ~ burn_p, disc_choli, mean) # en promedio es mejor la p grande




# Cómo evaluar el efecto de achicar el paisaje?? 
# quizás agarrar los fuegos grandes (p = 0.28) y apagarle todos los pixeles,
# imitando un paisaje más chico. luego volver a calcular los overlaps.

# Esto mismo se puede hacer sin fuegos, únicamente con el paisaje, usando distintos
# tamaños y quemando todo lo quemable.

# Problema hipotético: si el paisaje es chico, muchos fuegos pueden quedarse sin
# paisaje. El tema es que esos fuegos tienen que tener tan mal overlap como 
# los fuegos malos que sean pequeños. El problema de paisajes chicos es que 
# pueden suavizar el mal ajuste de los parámetros que hacen que todo se vaya a 
# la mierda.

# quizás estoy teniendo un problema visual. Con fuegos grandes es muy evidente 
# el mal overlap, pero con fuegos peques no es tan notable. Pero los números no 
# mienten... debería creerle más a los números.


# landscape size thoughts -------------------------------------------------

# overlap = intersection / union 

# Some thoughts imagining nested fires...

# With the smallest observed fire (10 ha), the smallest simulated fire (1 pixel)
# would take this overlap:

# overlap = 1 pixel / 10 ha pixels

# 1 ha = 10000 m2
# 1 pixel = 900 m2 (30 * 30) = 900 / 10000 = 0.09 ha

# in pixels, 10 ha is 10 * 10000 / 900 = 111.11111 pixels

# overlap_smallest_extreme = 0.09 / 10            = 0.009  # from ha
#                          = 1    /  111.11111111 = 0.009  # from pixels

# so, with the smallest fire, our resolution (30 m) does not allow for simulated
# fires on the small side to lower the overlap below 0.009.

# imagining a nested fire that burns the whole landscape, which would the fire
# size be?
# overlap_largest_extreme = 10 / X = 0.009
#                       X = 10 / 0.009 = 1111.111 ha

# So, to be symmetric, the relation of the largest to observed fire should be
# around
# 1 / 111.11 = 0.009               

# Cholila size in pixels:       sum(land[, "burned"])   = 317956
# Cholila landscape (burnable): sum(land[, "burnable"]) = 7122089
# Cholila landscape (all):      nrow(land)              = 10434109
                       
# Ratio cholila / burnable = 317956 / 7122089  =     0.04464364  >> 0.009
# Ratio cholila / all      = 317956 / 10434109 =     0.03047275  >> 0.009

# What if the landscape would be reduced to 2/3 ?
# new_land = 10434109 * (2/3) = 6956073

# Ratio cholila / all      = 317956 / 6956073  =     0.04570912  >> 0.009
# Expecting the same ratio o burnable:

# new_land_burnable = 10434109 * (2/3) * (7122089 / 10434109) = 4748059

# Ratio cholila / burnable = 317956 / 4748059  =     0.06696547  >> 0.009 

# and a disjoint fire with the same size?
# imagine
# (317956 * 0.01) / (317956 * 2 - (317956 * 0.01)) = 0.005025126

# So, a fire that burns all the landscape would be better than a disjoint one.
# Perhaps that makes sense. I've seen that the model with p = 0.28 generates
# small fires and fires that burn the whole landscape, so the variation is 
# expected.

# What follows? Make some benchmarks comparing large fires (p = 1 for cholila)
# and fires with constrained parameters for wind, aspect, elevation and slope.

# Theoretically, if only positive parameters are allowed, the landscape should
# constrain more the spread, making fires smaller.

# Then, perhaps benchmarking another landscape could be useful, as it could
# change dramatically the timing, and perhaps landscape size is not a so big
# problem (maybe just clipping cholila landscape?).

# Also, read about Gaussian process fitting. How many particles are used by wave?

