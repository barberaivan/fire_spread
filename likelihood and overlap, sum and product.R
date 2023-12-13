library(tidyverse)
normalize <- function(x) x / sum(x)

x <- seq(-5, 5, length.out = 1000)

d2 <- (dnorm(x) * 5) %>% normalize
dquad <- (dnorm(x) ^ 5) %>% normalize

plot(dquad ~ x, type = "l")
lines(d2 ~ x, col = "red")


setwd("/home/ivan/Insync/Fire spread modelling/fire_spread/files/pseudolikelihood_estimation")
ff <- list.files()
ff2 <- ff[grep("gp-fitted", ff, invert = T)]

tabs <- lapply(ff2, function(x) readRDS(x))
lapply(tabs, dim)

aa <- abind::abind(tabs, along = 1)[, "overlap", ]
dim(aa)
logsum <- apply(log(aa), 1, sum)
range(logsum)
logsum_shift <- logsum - max(logsum) # shift in log to get numbers 
# further away from zero

range(logsum_shift)

ss <- apply(aa, 1, mean)


ppn <- normalize(exp(logsum_shift)) * 100
ssn <- normalize(ss)  * 100

summary(ppn); summary(ssn)

plot(ppn[order(ppn, decreasing = T)])
plot(ssn[order(ssn, decreasing = T)])

plot(logsum_shift)
plot(ppn)
plot(aa[which.max(ppn), ])
plot(aa[6, ])


# try to fit the gp at the scale of logsum or at the scale of 
# product. that could make a huge difference. 