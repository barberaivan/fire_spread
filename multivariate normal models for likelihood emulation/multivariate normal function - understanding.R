# The multivariate normal density has this shape:
# exp( - (X - MU) %*% P %*% (X - MU) ) 

# Taking its log, we have
# - (X - MU) %*% P %*% (X - MU)

# So, it's just a quadratic shape.
# We must set initial values so that whithin the range of the predictor
# (X) the density goes from 0 to log(lower_similarity ~ 0.15) = -1.89712

log_simil_lower <- log(0.15)
lower_x <- -3
upper_x <- 3

s <- 2 # 0.1 o 0.2 a 20 parece bastante razonable, considerando a x en [-3, 3]
a <- 1 / s ^ 2 
curve(-a * x ^ 2, from = lower_x, to = upper_x, ylim = c(log_simil_lower - 1, 0.1),
      n = 1000)
abline(h = log_simil_lower, col = 2)





# skew normal?
# dnorm(x) * pnorm(slnat * x)
# log(dnorm(x)) + log(pnorm(alpha * x))
alpha = 10
curve(2 * dnorm(x) * 1/2 * pnorm(x * alpha), from = -3, to = 3)
abline(v = 0)

# log scale
alpha = -0
curve(dnorm(x, log = TRUE) + pnorm(x * alpha, log.p = TRUE), 
      from = -3, to = 3,
      ylim = c(-3, 0), n = 1000)
abline(v = 0, h = log(0.1), lty = 2, col = 2)


# DEFINITIVAMENTE ESTUDIAR LA POSTERIOR DE ESTOS MODELOS ME VA A AYUDAR. 

# PROBAR LA ESTIMACIÓN DE UNA NORMAL UNIVARIADA EN STAN PARA VER PROBLEMAS 
# DE CORR Y ESO. 

# LUEGO PASAR AL CASO BIVARIADO Y VER QUÉ ONDA, SI HAY PROBLEMAS. 

# LUEGO PASAR AL CASO UNIVARIADO SKEW.

# LUEGO, MULTIVARIADO SKEW.

# Betancourt tiene un blog sobre GP que está terrible:
# Ahí encuentra que la likelihood hace cualquiera y que hay que penalizarla
# para que de resultados sensatos. Probablemente lo mismo esté pasando acá.
# Incluso si quiero usar estimación puntual, necesito regularizar la likelihood,
# o sea, estimar el posterior mode.
# https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#23_Fitting_A_General_Gaussian_Process_Posterior


# Funciones de stan para calcular ahí dentro las densities
# https://mc-stan.org/docs/functions-reference/multivariate-normal-distribution.html
# https://mc-stan.org/docs/functions-reference/normal-distribution.html