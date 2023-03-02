# starting with spread function

# para llamar a pak (renv no me deja)
library(bench)
library(terra)
library(tidyverse)

# spread function 
spread_prob <- function(focal, targets) {
  
  # focal = burning_cells[b]
  # targets = neighbours
  
  data_targets <- values(r_predictors)[targets, , drop = FALSE]
  data_focal <- values(r_predictors)[focal, ]
  
  # recompute wind and elev columns (will be slope and wind effects)
  data_targets[, "elev"] <- sin(atan((data_targets[, "elev"] - data_focal["elev"]) / d[cols_use]))
  data_targets[, "wind"] <- cos(angles[cols_use] - data_focal["wind"])
  
  # careful with upper limit
  probs <- plogis(coefs[1] + data_targets %*% coefs[2:length(coefs)]) * 0.5
  return(probs)
}


# model coefficients (includes intercept)
coefs <- c(0.5, -1.5, 0.5, 0.5, 0.5, 0.5, 0.5)

# test wind effect
coefs <- c(0.5, 0, 0, 0, 3, 0, 0)

# test slope effect
coefs <- c(0.5, 0, 10, 0, 0, 0, 0)

names(coefs) <- c("Intercept", "veg", "elev", "aspect", "wind", "pp", "temp")

# predictors raster

size <- 50

n_rows <- size
n_cols <- size

r_veg <- rast(ncol = n_cols, nrow = n_rows,
              xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
r_elev <- r_veg
r_aspect <- r_veg
r_wind <- r_veg
r_pp <- r_veg
r_temp <- r_veg

values(r_veg) <- rbinom(ncell(r_veg), prob = 0.5, size = 1)
# to test elevation:
elevs <- matrix(NA, nrow(r_veg), ncol(r_veg))
elevs[1, ] <- 2000
for(i in 2:nrow(elevs)) elevs[i, ] <- elevs[i-1, ] - 2000 / nrow(r_veg)

values(r_elev) <- as.numeric(t(elevs))#rnorm(ncell(r_veg), 1000, 100) %>% abs
values(r_aspect) <- cos(runif(ncell(r_veg), 0, 2 * pi) - 315 * pi / 180)
values(r_wind) <- 0#pi#runif(ncell(r_veg), 0, 2 * pi)
values(r_pp) <- runif(ncell(r_veg), 0, 1)
values(r_temp) <- runif(ncell(r_veg), 0, 1)


r_predictors <- c(r_veg, r_elev, r_aspect, r_wind, r_pp, r_temp)
names(r_predictors) <- names(coefs)[-1]

# values(r_predictors)[10, ]

# vector of distances and angles between a cell and its neighbours
# (used to compute slope and wind effects)

d <- rep(res(r_predictors)[1], 8) # sides
d[c(1, 3, 6, 8)] <- res(r_predictors)[1] * sqrt(2)

angles_matrix <- matrix(
  c(315, 0, 45,
    270, NA, 90,
    225, 180, 135), 
  3, 3, byrow = TRUE) * pi / 180
angles <- as.numeric(angles_matrix %>% t) %>% na.omit

# create example fire raster

r <- r_veg
# fill raster
# 0: unburned
# 1: burning
# 2: burned
# 3: not flammable
values(r) <- 0 #rpois(ncell(r), lambda = 1)

# ignition point:
ig_location <- sample(1:ncell(r), size = 1, replace = FALSE) #floor(ncell(r) / 2 - 2) # cell number
#ig_location <- cellFromRowCol(r, nrow(r) / 2, ncol(r) / 2)
#plot(r_elev)
# set as burning
values(r)[ig_location] <- 1


# initialize burning_cells cell id
burning_cells <- which(values(r) == 1) # cell number id
j = 1
plot(r, main = paste("step = ", j), col = c("green", "red"))

# spread
while(length(burning_cells) > 0) {
  
  print(paste("cycle ", j, sep = ""))
  
  # update: burning to burned
  values(r)[burning_cells] <- 2

  # get cell id from neighbours
  neighbours_matrix <- adjacent(r, cells = burning_cells, "queen") 
  #colnames(neighbours_matrix) <- 1:8
  
  # spread from burning pixels
  for(b in 1:length(burning_cells)) {
    print(paste("burning cell ", b, sep = ""))
    
    #b = 1
    # get neighbours available to burn (cell ids), while getting rid of 
    # non-existent neighbours
    
    filter <- !is.na(values(r)[neighbours_matrix[b, ]]) &
              values(r)[neighbours_matrix[b, ]] == 0
    cols_use <- which(filter)
    neighbours <- neighbours_matrix[b, cols_use]
    
    # simulate spread
    if(length(neighbours) > 0) {
      
      ##### compute probs
      
      # focal = burning_cells[b]
      # targets = neighbours
      # 
      # data_targets <- values(r_predictors)[targets, , drop = FALSE]
      # data_focal <- values(r_predictors)[focal, ]
      # 
      # # recompute wind and elev columns (will be slope and wind effects)
      # data_targets[, "elev"] <- sin(atan((data_targets[, "elev"] - data_focal["elev"]) / d[cols_use]))
      # data_targets[, "wind"] <- cos(angles[cols_use] - data_focal["wind"])
      # 
      # # careful with upper limit
      # probs <- plogis(coefs[1] + data_targets %*% coefs[2:length(coefs)]) * 0.5
      
      #####
      
      values(r)[neighbours] <- rbinom(
        n = length(neighbours), size = 1,
        #prob = probs
        prob = spread_prob(focal = burning_cells[b], targets = neighbours)
      )
    }
    
  }
  
  j <- j + 1
  plot(r, main = paste("step = ", j), col = c("green", "red", "black"))
  # update burning_cells
  burning_cells <- which(values(r) == 1)
}
#plot(r, main = paste("step = ", j), col = c("green", "red", "black"))

# hay que comparar con la forma bruta que evalúa todos los píxeles aunque ya se
# hayan quemado, que no es recursiva.

# También hay que probar escribir con rcpp todo lo que se pueda, a ver si acelera.
# Incluso se podría escribir todo en rcpp si pasamos todos los datos a algo matricial.

plot(r)
mmm <- as.matrix(r)
str(mmm)



# subsetting and getting coordinates

rowFromCell(r, 150)
## [1] 3
colFromCell(r, 150)
## [1] 28
cellFromRowCol(r,5,5)
## [1] 149
xyFromCell(r, 100)
##       x  y
## [1,] 95 65
cellFromXY(r, cbind(0,0))
## [1] 343
colFromX(r, 0)
## [1] 19
rowFromY(r, 0)
## [1] 10


# angles and radians
315 * pi / 180
135 * pi / 180

(135-315) * pi / 180

cos(-pi)
cos(pi*1.5)
cos(0)

# angulos around a cell
seq(0, 360, by = 360 / 8)[1:8]
seq(0, 360, by = 360 / 8)[1:8] * pi / 180
ang_neigh <- seq(0, 360, by = 360 / 8)[1:8] * pi / 180
cos(runif(1, 0, 2 * pi) - ang_neigh)

sin(45 * pi / 180)
sin(90 * pi / 180)

pi / 2 # = 1.57
atan(10000 / 10) # devuelve el ángulo en radianes

delta_alt <- 30
dist_cell <- 30
sin(atan(delta_alt / dist_cell))

atan(delta_alt / dist_cell) * 180 / pi # perfect, 45 °

sqrt(30 ^ 2 + 30 ^ 2); 30 * sqrt(2)



# Ideas -------------------------------------------------------------------

# terra está escrito en C++, así que no tendría mucho sentido pasar sus 
# funciones a Rcpp. 
# Así, sería sensato seguir trabajando de forma similar, pasando todo lo que 
# pueda a rcpp pero sin pasar las funciones de terra, como adjacent.

# Temas de memoria:
# raster de tamaño around cholila:
r_large <- rast(ncol = 1900, nrow = 1900,
                xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
values(r_large) <- 0
r_large
object.size(r_large)
adj_large <- adjacent(r_large, 1:ncell(r_large), direction = "queen")
# este pesa 462 Mb, así que no sería bueno precalcular todas las adyacencias.
object.size(adj_large)
object.size(cbind(values(r_large), adj_large))
class(adj_large)

# pruebo cuánto pesa una large image
library(terra)
dnbr <- rast("/home/ivan/Insync/Burned area mapping/Long mapping approach/burned_area_context_100m_connected_components.tif")
inMemory(dnbr) # NOT
#dnbr_val <- values(dnbr) # esto hace abortar la misión
values(dnbr)[100, 1] # también explota la memoria
# este raster es terriblemente extremo, pero igual tendremos alrededor de 7 
# capas para cada fuego, y sumando fuegos quizás lleguemos al tope de RAM. 

# Sobre memoria

# Mejor primero generar los verdaderos rasters para cada fuego y ahí ver 
# si se puede sostener todo en RAM. Si es posible, pasarse a Rcpp.
# O comparar el uso de memoria de los rasters y las matrices.
# Es raro, parece que la raster en sí ocupa poco lugar, pero al pedir sus valores
# usa mucho. 

# Spread function options:

# Cuando evalúo los targets uno por uno, de forma recursiva y escribiendo el raster
# apenas calculo la igninción (o su ausencia), tengo que acceder al raster muchas 
# veces. Pero si voy guardando eso y accedo al raster una única vez, quizás ahorre 
# tiempo.

# prueba pequeña:
r <- rast(ncol = 10, nrow = 30,
          xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
values(r) <- 0

fonce <- function(x) {
  values(r)[1:300] <<- rnorm(300)
  return("done")
}

floop <- function(x) {
  for(i in 1:300) values(r)[i] <<- rnorm(1)
  return("done")
}

bench::mark(
  once = fonce(1),
  loop = floop(1),
  iterations = 5
)

fonce()

system.time(
  values(r)[1:300] <- rnorm(300)
)
system.time(
  for(i in 1:300) values(r)[i] <- rnorm(1)
)

# definitivamente, acceder solo una vez es mucho más rápido.

# Comparar usando solo R functions.
# Eventualmente, todo lo que se hace dentro de un update puede estar dentro de 
# Rcpp.
# Antes y después de eso hay que subsetear valores del raster, y a eso lo hago 
# con terra.

# update once by iter -------------------------------------------------------

# vector of distances and angles between a cell and its neighbours
# (used to compute slope and wind effects)

d <- rep(res(r_predictors)[1], 8) # sides
d[c(1, 3, 6, 8)] <- res(r_predictors)[1] * sqrt(2) # diagonal distances

angles_matrix <- matrix(
  c(315, 0, 45,
    270, NA, 90,
    225, 180, 135), 
  3, 3, byrow = TRUE) * pi / 180
angles <- as.numeric(angles_matrix %>% t) %>% na.omit

# model coefficients
coefs <- c(3, 0, 0, 0, 0, 0, 0)
names(coefs) <- c("Intercept", "veg", "elev", "aspect", "wind", "pp", "temp")

# Create rasters
size <- 100
n_rows <- size
n_cols <- size

# fire raster
r <- rast(ncol = size, nrow = size,
          xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)

# predictors
r_veg <- r
r_elev <- r
r_aspect <- r
r_wind <- r
r_pp <- r
r_temp <- r

values(r_veg) <- rbinom(ncell(r), prob = 0.5, size = 1)
values(r_elev) <- rnorm(ncell(r), 1000, 100) %>% abs
values(r_aspect) <- cos(runif(ncell(r), 0, 2 * pi) - 315 * pi / 180)
values(r_wind) <- 0#pi#runif(ncell(r_veg), 0, 2 * pi)
values(r_pp) <- runif(ncell(r), 0, 1)
values(r_temp) <- runif(ncell(r), 0, 1)

r_predictors <- c(r_veg, r_elev, r_aspect, r_wind, r_pp, r_temp)
names(r_predictors) <- names(coefs)[-1]

# set fire values
values(r) <- 0

# ignition point:
burning_cells <- 2#floor(ncell(r) / 2) - ncol(r) / 2#sample(1:ncell(r), size = 1, replace = FALSE) #floor(ncell(r) / 2 - 2) # cell number
burning_cells <- burning_cells:(burning_cells + 1)

# spread
j = 1
while(length(burning_cells) > 0) {
  
  # visualize
  values(r)[burning_cells] <- 1 # only needed to plot
  plot(r, main = paste("step = ", j), col = c("green", "red", "black"))
  # update: burning to burned (done here only to plot)
  values(r)[burning_cells] <- 2
  
  print(paste("cycle ", j, sep = ""))

  # get cell id from neighbours
  neighbours_matrix <- adjacent(r, cells = burning_cells, "queen") 
  #neighbours_matrix <- adjacent(r, cells = burning_cells, "queen", pairs = TRUE) 
  
  neigh_val <- values(r_predictors)[neighbours_matrix, ]
  # a neigh_val le agrego 3 columnas: 
  # - posición (1:8), indicando qué vecino es
  # - focal (burning_cells), indicando de qué celda prendida viene (1:length(burning_cells)) y
  # - id, la celda a la que corresponde.
  
  extrad <- matrix(c(rep(1:8, each = length(burning_cells)),
                     rep(1:length(burning_cells), 8),
                     as.numeric(neighbours_matrix),
                     rep(0, length(burning_cells) * 8)),
                   ncol = 4)
  colnames(extrad) <- c("position", "focal", "id", "burn")
  
  # identify burnable neighbours
  neigh_use <- !is.na(extrad[, "id"]) & values(r)[extrad[, "id"]] == 0
  ids_use <- (1:nrow(extrad))[neigh_use]
  
  # get values from burning_cells
  focal_val <- values(r_predictors)[burning_cells, , drop = F]
  
  for(i in ids_use) {
    #i = 7
    focal = extrad[i, "focal"]
    
    # recompute wind and elev columns (will be slope and wind effects)
    neigh_val[i, "elev"] <- sin(atan((neigh_val[i, "elev"] - focal_val[focal, "elev"]) / d[extrad[i, "position"]]))
    neigh_val[i, "wind"] <- cos(angles[extrad[i, "position"]] - focal_val[focal, "wind"])

    # careful with upper limit
    extrad[i, "burn"] <- rbinom(n = 1, size = 1, 
                                prob = plogis(coefs[1] + neigh_val[i, , drop = F] %*% coefs[2:length(coefs)]) * 0.4)
  }
  
  # get new burning_cells
  burning_cells <- extrad[extrad[, "burn"] == 1, "id"] %>% unique
  
  # advance step
  j <- j + 1

}



# comparando el object_size de raster y array

object.size(r_predictors)
r_pred_val <- values(r_predictors)
object.size(r_pred_val)
r_predictors

object.size(r_predictors@ptr$valueType)
names(r_predictors@ptr$.pointer)

# no puedo comparar porque el raster solo tiene un pointer a un objeto de c++

class(r_predictors)

