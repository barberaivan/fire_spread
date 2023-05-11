# The aspect importance (in [0, 1]) should be an increasing function of slope, 
# passing through the origin. A piece-wise linear or a monomolecular function 
# could be used, but the functions require parameterization.

spat_data <- read.csv("/home/ivan/Insync/patagonian_fires/burned_area_patterns/data_spatial_variables.csv")
spat_data <- spat_data[spat_data$class == "Available", ]

slope <- spat_data$slope
quantile(slope, probs = seq(0, 1, by = 0.1))["80%"]
secdf <- ecdf(slope)
plot(secdf)
secdf(28)
secdf(15)

curve(1 - exp(-8 * ((x - 5) * pi/180)), from = 0, to = 80)
abline(v = c(5, 28), h = 1)

curve((1/28) * x, from = 0, to = 28)
abline(v = c(5, 28), h = 1)


curve((1/28) * x, from = 0, to = 28)
abline(v = c(5, 28), h = 1)

curve(plogis(-5 + 0.5 * x), from = 0, to = 45, n = 1000)
abline(v = c(5, 28), h = 0)

aspect_importance <- function(x, a = -5, b = 0.5) {
  plogis(a + b * x)
}

imp <- aspect_importance(slope)

# compute linear aspect (northing)
aspect <- cos(spat_data$aspect * pi / 180)

a <- -5
b <- 0.35
weight <- aspect_importance(slope, a, b)

par(mfrow = c(2, 3))
hist(slope, main = "slope distribution", xlab = "slope (°)")
abline(v = c(5, 28), lty = 2, col = 2)
curve(plogis(-5 + 0.35 * x), from = 0, to = 45, n = 1000,
      ylab = "weight", xlab = "slope (°)",
      main = "northing importance function",
      ylim = c(0, 1))
abline(v = c(5, 28), h = c(0, 1), lty = 2, col = 2)
hist(weight, main = "weight distribution")
# hist(slope)
hist(aspect, main = "northing distribution",
     xlab = "northing")
hist(aspect * weight, main = "weighted northing distribution",
     xlab = "northing * weight")
par(mfrow = c(1, 1))


par(mfrow = c(2, 3))
plot(density(slope, from = 0), 
     main = "slope distribution", xlab = "slope (°)")
abline(v = c(5, 28), lty = 2, col = 2)
curve(plogis(-5 + 0.35 * x), from = 0, to = 45, n = 1000,
      ylab = "weight", xlab = "slope (°)",
      main = "northing importance function",
      ylim = c(0, 1))
abline(v = c(5, 28), h = c(0, 1), lty = 2, col = 2)
plot(density(weight, from = 0, to = 1), 
     main = "weight distribution")
# plot(density(slope)
plot(density(aspect, from = -1, to = 1), 
     main = "northing distribution",
     xlab = "northing")
plot(density(aspect * weight, from = -1, to = 1, adjust = 0.5), 
     main = "weighted northing distribution",
     xlab = "northing")
par(mfrow = c(1, 1))

