# After separating the fires with many parts, run the following:

library(terra)
library(tidyverse)

# get path for data
local_dir <- normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
dir_split <- strsplit(local_dir, .Platform$file.sep)[[1]]
# replace the "fire_spread" directory by "data"
dir_split[length(dir_split)] <- "fire_spread_data"
data_path <- paste(dir_split, collapse = .Platform$file.sep)

lands <- readRDS(file.path(data_path, "landscapes_ig-known_non-steppe.rds"))

mean_sd_rast <- function(arr) {
  fwi <- as.numeric(arr[, , "fwi"])
  m <- mean(fwi)
  s <- sd(fwi)
  a <- sum(as.numeric(arr[, , "burned"]))
  return(c("mean" = m, "sd" = s, "area_pix" = a))
}

summs <- sapply(lands, mean_sd_rast) %>% t %>% as.data.frame

years <- sapply(names(lands), function(n) {
  strsplit(n, "_")[[1]][1] 
}) %>% unname
years[years == "CoVentana"] <- "1999"
years[years == "CoCampana"] <- "2002"

summs$year <- as.numeric(years)
summs$area_ha <- summs$area_pix * 0.09 # pixel is 30 * 30 m


# Plots of Fire size as a function of FWI

theme_set(theme_bw())

ggplot(summs, aes(mean, area_ha)) + 
  geom_smooth(method = "glm",
              method.args = list(family = Gamma(link = "log"))) +
  # geom_smooth() +
  geom_point(size = 3, alpha = 0.7) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Fire area (ha)") + 
  xlab("FWI standardized (landscape mean)")

ggsave("files/fwi_area.png",
       width = 10, height = 8, units = "cm")

ggplot(summs, aes(mean, area_ha)) + 
  geom_smooth(method = "lm") + 
  geom_point(size = 3, alpha = 0.7) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Fire area (ha, log10)") + 
  xlab("FWI standardized (landscape mean)") +
  scale_y_continuous(trans = "log10")

ggsave("files/fwi_area_log10.png",
       width = 10, height = 8, units = "cm")


# FWI intra- and inter-years variation -----------------------------------

ydata <- aggregate(cbind(mean, sd) ~ year, summs, mean)
hist(ydata$sd, xlim = c(0, 1.5)); abline(v = sd(ydata$mean))

hist(summs$sd, xlim = c(0, 2),
     xlab = "sd(FWI) whithin landscapes", 
     main = NULL)
abline(v = sd(ydata$mean), col = "red", lty = 2, lwd = 1.5)
abline(v = sd(summs$mean), col = "blue", lty = 2, lwd = 1.5)
text(0.6, 35, "sd(FWI)\nbetween years", col = "red")
text(1.6, 35, "sd(FWI)\nbetween fires", col = "blue")


plot(sd ~ area_ha, summs)
