# Pidiendo 90 m de res para cholila anduvo bien.
# Pidiendo 30 m, crashea.

library(terra)

angles <- rast("/home/ivan/Insync/Fire spread modelling/data/focal fires data/data_cholila_elevation_284_10_90m_ang.asc")
angles
res(angles)
plot(angles)

# Interpolation to 30 m (we will use just a projection instead)

# bring elevation data (30 m raster)
elev <- rast("/home/ivan/Insync/Fire spread modelling/data/focal fires data/data_cholila_elevation.tif")

angles30 <- project(angles, elev, method = "cubicspline")
plot(angles30)

size(angles); size(angles30); size(elev) # OK
res(angles30)
