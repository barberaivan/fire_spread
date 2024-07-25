library(tidyverse)
library(mgcv)
library(viridis)

# Square function ---------------------------------------------------------

N <- 2000
x <- runif(N)
z <- runif(N)
y <- as.numeric(x > 0.5 & z < 0.8)
df <- data.frame(x = x, y = y, z = z)

m1 <- gam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr"),
          method = "REML", family = "binomial")

pd <- expand.grid(x = seq(0, 1, length.out = 100),
                  z = seq(0, 1, length.out = 100))
pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.5,
             color = viridis(1, begin = 0.5)) +
  scale_fill_viridis(option = "B")


# Square at centre function -----------------------------------------------

N <- 2000
x <- runif(N)
z <- runif(N)
y <- as.numeric(x > 0.5 & x < 0.8 & z < 0.8 & z > 0.5)
df <- data.frame(x = x, y = y, z = z)

m1 <- gam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr"),
          method = "REML", family = "binomial")

pd <- expand.grid(x = seq(0, 1, length.out = 100),
                  z = seq(0, 1, length.out = 100))
pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")

# Doesnt need interaction


# Triangle ---------------------------------------------------------------

N <- 2000
x <- runif(N)
z <- runif(N)
y <- as.numeric(x > 0.5 * z + 0.3)

df <- data.frame(x = x, y = y, z = z)


pd <- expand.grid(x = seq(0, 1, length.out = 100),
                  z = seq(0, 1, length.out = 100))

m1 <- bam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr"),
          method = "REML", family = "binomial")

pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")

# Doesnt need interaction

# Fringe -----------------------------------------------------------------

N <- 2000
x <- runif(N)
z <- runif(N)
y <- as.numeric(x > z - 0.2 & x < z + 0.1)
df <- data.frame(x = x, y = y, z = z)

pd <- expand.grid(x = seq(0, 1, length.out = 100),
                  z = seq(0, 1, length.out = 100))

m1 <- bam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr"),
          method = "REML", family = "binomial")

pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")


m1 <- bam(y ~ te(x, z, k = 10, bs = "cr"),
          method = "REML", family = "binomial")
pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")


m1 <- bam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr") +
              ti(x, z, k = 10, bs = "cr"),
          method = "REML", family = "binomial")
pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")
summary(m1)

# The fringe needs interaction; squares or triangles, don't
# (as long as the triangle is in a corner)


# non-corner Triangle ---------------------------------------------------------------

N <- 2000
x <- runif(N)
z <- runif(N)
y <- as.numeric(x > z & z > 0.5 & x < 0.8 & x > 0.5)
df <- data.frame(x = x, y = y, z = z)

pd <- expand.grid(x = seq(0, 1, length.out = 100),
                  z = seq(0, 1, length.out = 100))

m1 <- bam(y ~ s(x, k = 15, bs = "cr") + s(z, k = 15, bs = "cr"),
          method = "REML", family = "binomial")

m1 <- bam(y ~ te(x, z, k = 10, bs = "cr"),
          method = "REML", family = "binomial")

pd$p <- predict(m1, pd, type = "response")

ggplot() +
  geom_tile(aes(x, z, fill = p), data = pd) +
  geom_point(aes(x, z), data = df[df$y == 1, ], alpha = 0.8,
             color = viridis(1, begin = 0.7)) +
  scale_fill_viridis(option = "B")


# When there are no fringes, the marginal basis work quite well. But triangles
# at the centre at better represented by interactions.