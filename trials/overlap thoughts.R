size <- 500
sims <- 1:(size * 6)
drel <- abs(sims - size) / size
plot(drel ~ sims)

ovsize <- 1000
unsize <- ovsize + 1:(ovsize * 100)

plot(ovsize / unsize ~ unsize)
plot((ovsize - unsize) / unsize ~ unsize)

plot((ovsize - unsize) ~ unsize)
plot((unsize - ovsize) / ovsize ~ unsize)

# overlap_delta
overlap_d <- (ovsize - unsize) / ovsize
overlap_q <- ovsize / unsize

plot(overlap_d ~ overlap_q, type = "l")
abline(v = 0.05, lty = 2)

# el overlap delta mira un efecto lineal. El cociente es no lineal.
# Entonces, overlaps_q de ~0.33 implican un delta = 2
approx(overlap_d, overlap_q, xout = -2)

approx(overlap_q, overlap_d, xout = 0.2)
# overlap de 0.2 implica un delta de 4.

# el overlap_q es ultrasensible en la escala de ~0.1 a 1
approx(overlap_q, overlap_d, xout = 0.1)
# o sea, diferencias de tamaño de 9 a 0. Pero varía a una escala muy pequeña
# en deltas > 9

# el paisaje suele estar limitado around delta =
approx(overlap_q, overlap_d, xout = 0.05) # -19

# si tenemos el overlap_q, el size y el size_diff, deberíamos poder llegar a
# overlap_delta, no?

overlap_q = common / (a + b - common)

a = fire_size
b = fire_size + size_diff

overlap_q * (a + b - common) = common
overlap_q * a + overlap_q * b - overlap_q * common = common
overlap_q * a + overlap_q * b = common + common * overlap_q
overlap_q * a + overlap_q * b = common * (1 + overlap_q)
common = (overlap_q * a + overlap_q * b) / (1 + overlap_q)

# más fácil sería
(overlap_q - 1) / overlap_q

## nooooo, pero así justamente most overlap_d irían como de -20 a -100, como
## si fueran un log, mucho peor que un logit. Justamente el overlap_q en escala
## original (no logit) es bastante insensible a las diferencias muy grandes,
## porque en la escala de las diff, lo que va de 20 a Inf, en el q va de
## 0.05 a 0.

