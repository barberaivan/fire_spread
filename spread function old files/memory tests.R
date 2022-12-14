library(terra)
library(tidyverse)

# start at 3.1 Gb
size <- 2000
r <- rast(ncol = size, nrow = size,
          xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
r
values(r) <- rnorm(size * size)
r
# llegamos a 3.2 Gb

m <- rnorm(size * size * 6)
format(object.size(m), units = "Mb")

m <- matrix(rnorm(size * size * 6), ncol = 6)
format(object.size(m), units = "Mb")

adj <- adjacent(r, cells = 1:ncell(r), directions = "queen")
format(object.size(adj), units = "Mb") # 488 Mb
class(adj)
str(adj)

m8 <- matrix(rnorm(size * size * 8), ncol = 8)
format(object.size(m8), units = "Mb") # 244 Mb, la mitad de lo que ocupa la matriz de adyacencia.

# Para que use menos ram debería escribir en c++ adjacent.

