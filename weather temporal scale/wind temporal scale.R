# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(terra)
library(lubridate)
library(viridis)
library(bayestestR)
library(circular) # density.circular, for wind direction map

# Functions ---------------------------------------------------------------

# mean, median and maximum
summ <- function(x) {
  return(
    c("mean" = mean(x),
      "mode" = map_estimate(x))
  )
}

summ3 <- function(x) {
  return(
    c("mean" = mean(x),
      "median" = median(x),
      "mode" = map_estimate(x))
  )
}

# cicular mean
# from
# https://stackoverflow.com/questions/32404222/circular-mean-in-r
mean_circular_raw <- function (x) { # takes angle in radians
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  return(circmean)
}

mean_circular_deg <- function (x) { # takes angle in degrees
  # x <- c(180, 160, 170)
  conv <- 2 * pi / 360 # degrees to radians factor
  mm_raw <- mean_circular_raw(conv * x) / conv
  mm <- (mm_raw + 360) %% 360
  return(mm)
}

summ_circ <- function(x, bw = 100) { # x in degrees
  # test
  # x <- rnorm(15, 270, 10); bw = 100

  # compute mode
  xc <- circular(x, type = "angles", units = "degrees", template = "geographic")
  dd <- density.circular(xc, bw = bw)

  # arrange angle scale
  den_x <- as.numeric(dd$x)
  den_x[den_x < 0] <- den_x[den_x < 0] + 360

  # mode
  x_mode <- den_x[which.max(as.numeric(dd$y))]

  res <- c("mean" = mean_circular_deg(x), "mode" = x_mode)
  return(res)
}

# wind direction equation from:
# https://stackoverflow.com/questions/21484558/how-to-calculate-wind-direction-from-u-and-v-wind-components-in-r
winddir_compute <- function(u, v) { # u = westwind (comes from west), v = southwind (comes from south)
  wind_abs = sqrt(u ^ 2 + v ^ 2)
  wind_dir_trig_to = atan2(u / wind_abs, v / wind_abs)
  wind_dir_trig_to_degrees = wind_dir_trig_to * 180 / pi

  # convert this wind vector to the meteorological convention of the direction
  # the wind is coming from:
  wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180

  return(wind_dir_trig_from_degrees)
}



# Data --------------------------------------------------------------------

fires_vec <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
fires <- as.data.frame(fires_vec)
# remove fires with uncertain year
uncertain <- which(fires$obs == "uncertain_year")
fires <- fires[-uncertain, ]

wind_data <- read.csv("data/climatic_data_by_fire_wind.csv")
# windspeed u and v components are in m/s

# Turn dates into dates
dd <- c("date_start", "date_end")
for(d in dd) wind_data[, d] <- as.Date(wind_data[, d])

dd <- c("date", "date_l", "date_u")
for(d in dd) fires[, d] <- as.Date(fires[, d])

# Correct date_l and date_u for Landsat-based images, so they are included in the
# focal period. For landsat, date_l and date_u are the dates between which the
# fire occurred, but not included. For modis, they are inclusive. Turn all into
# inclusive. Do this only for date_diff larger than 2
ll <- (fires$datesrc %in% c("landsat", "Mermoz")) &
  (fires$date_u - fires$date_l > 2)
fires$date_l[ll] <- fires$date_l[ll] + 1
fires$date_u[ll] <- fires$date_u[ll] - 1

# date uncertainty
fires$datediff <- fires$date_u - fires$date_l
rows_fit <- fires$datediff < 20

# fires without much date uncertainty
fires_sub <- fires[rows_fit, ]

# compare middate modis with thomas dates
fires$middate_sat <- fires$date_l + floor((fires$date_u - fires$date_l) / 2)

delta_sat_th <- fires$date - fires$middate_sat
# hist(as.numeric(delta_sat_th[delta_sat_th < 100]), breaks = 30)
# hist(as.numeric(delta_sat_th[grep("modis", fires$datesrc)]), breaks = 30)
# plot(table(abs(as.numeric(delta_sat_th[grep("modis", fires$datesrc)]))))
# plot(table(abs(as.numeric(delta_sat_th[abs(delta_sat_th) < 100]))))

# compute wind direction and speed
wind <- data.frame(
  fire_id = wind_data$fire_id,
  date_start = wind_data$date_start,
  date_end = wind_data$date_end
)
wind$westwind <- wind_data[, grep("west", names(wind_data))] %>% as.matrix
wind$southwind <- wind_data[, grep("south", names(wind_data))] %>% as.matrix
wind$speed <- sqrt(wind$westwind ^ 2 + wind$southwind ^ 2) * 3.6 # scale to km/h
wind$direction <- winddir_compute(wind$westwind, wind$southwind)

# plot(density(as.numeric(wind$speed))) # OK
# wcirc <- circular(wind$direction[1, 1:14], type = "angles", units = "degrees",
#                   template = "geographic")
# d_dir <- density.circular(wcirc, bw = 100)
# plot(d_dir)

# compute the wind summaries (summ) for fortnight and for focal date ------

# wind speed and direction for focal date and fortnight. For the latter,
# get two summaries: mean and map.

# get mean and mode aggregated over fortnite 
wind_agg_fort <- do.call("rbind", lapply(fires$fire_id, function(f) {
  # f = fires$fire_id[41]

  # get wind values
  row <- which(fires$fire_id == f)
  dtemp <- wind[wind$fire_id == f, , drop = F]

  #wind_vals <- dtemp[1, grep("wind", names(wind_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)

  # Choose middate here
  middate <- fires$date[row]
  # middate <- fires$middate_sat[row]

  # get week to get fortnight
  dd <- data.frame(direction = as.numeric(dtemp$direction),
                   speed = as.numeric(dtemp$speed),
                   date = date_seq_av,
                   week = week(date_seq_av),
                   year = year(date_seq_av))

  dd <- dd[order(dd$date), ]

  # if there is a year change in the period, restart weeks backwards from
  # december 31st. This is because the 7-day counting may create a last week
  # with too few days.
  if(max(diff(dd$year)) == 1) {
    restart <- dd$year == min(dd$year)
    weeks_n <- ceiling(sum(restart) / 7)
    week_tmp <- rep(seq(0, -(weeks_n - 1)), each = 7)[1:sum(restart)]
    week_use <- rev(week_tmp)
    dd$week[restart] <- week_use
  }

  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

  dd_mid <- dd[dd$date == middate, , drop = F]

  # obtain values
  speed_focal <- dd_mid$speed
  dir_focal <- dd_mid$direction

  speed_fort <- summ(dd$speed[dd$fortnight == dd_mid$fortnight])
  dir_fort <- summ_circ(dd$direction[dd$fortnight == dd_mid$fortnight])

  res_dir <- matrix(c(dir_focal, dir_fort), 1)
  res_speed <- matrix(c(speed_focal, speed_fort), 1)
  colnames(res_dir) <- colnames(res_speed) <- c("focal", "fort_mean", "fort_mode")

  # data frame to export result
  res <- fires[row, c("fire_id", "date")]
  res$direction <- res_dir
  res$speed <- res_speed

  return(res)
})) %>% as.data.frame()

# View(wind_agg_fort)

# GGally::ggpairs(wind_agg_fort$direction %>% as.data.frame())
# GGally::ggpairs(wind_agg_fort$speed %>% as.data.frame())

# save using only the mean and the map by fortnight


# compute the wind summaries (summ) for focal period ---------------------

# compute direction and speed based on the summarized u and v components.
# when date is quite certain (< 5 days range), use an extended period to 
# smooth easterly winds that are not consistent with most fires.

# The period to summarize is defined as follows:

# Date defined with records:
#   here date is usually the lowest bound, so take the timespan from 
#   date to date_u.
#   If date_u - date > 15 days, take
#   date + 6 days (one week).
#   If 0 <= date_u - date < 5, take
#   date + 6 days.
#   If date_u - date < 0, take
#   [date - 3, date + 3]

# Date not defined by records (modis or landsat):
#   If date_u - date_l >= 5, average between date_u and date_l.
#   If 0 <= date_u - date_l < 5 average between middate - 3 and middate + 3

records <- 1:nrow(fires) %in% grep("record", fires$datesrc)

wind_agg <- do.call("rbind", lapply(fires$fire_id, function(f) {
  # f = fires$fire_id[41]
  
  # get wind values
  row <- which(fires$fire_id == f)
  dtemp <- wind[wind$fire_id == f, , drop = F]
  
  # define available  dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  
  # define date range to summarize
  if(records[row]) {
    datediff <- fires$date_u[row] - fires$date[row]
    lower <- fires$date[row]
    
    if(datediff > 15) upper <- fires$date[row] + days(6)
    if(datediff <= 15 & datediff >= 5) upper <- fires$date_u[row]
    if(datediff < 5 & datediff >= 0) upper <- fires$date[row] + days(6)
    if(datediff < 0) {
      lower <- fires$date[row] - 3
      upper <- fires$date[row] + 3
    }
  } else { # not record-based dates
    datediff <- fires$date_u[row] - fires$date_l[row]
    if(datediff >= 5) {
      lower <- fires$date_l[row]
      upper <- fires$date_u[row]
    } 
    if(datediff < 5 & datediff >= 0) {
      middate <- fires$date_l[row] + days(ceiling(datediff / 2))
      lower <- middate - days(3) 
      upper <- middate + days(3) 
    }
    if(datediff < 0) {
      lower <- fires$date[row] - 3
      upper <- fires$date[row] + 3
    }
  }
  
  date_seq_summ <- seq(lower, upper, 1)
  
  # get wind data in those days
  ids <- which(date_seq_av %in% date_seq_summ)
  ww <- wind$westwind[row, ids] %>% unname
  sw <- wind$southwind[row, ids] %>% unname
  
  uv_summ <- data.frame(
    "west" = summ3(ww),
    "south" = summ3(sw)
  )
  
  winddir <- winddir_compute(uv_summ$west, uv_summ$south)
  windspeed <- sqrt(uv_summ$west ^ 2 + uv_summ$south ^ 2) * 3.6 # scale to km/h
  
  res <- matrix(c(winddir, windspeed), 1)
  colnames(res) <- paste(rep(c("direction", "speed"), each = 3),
                         rep(rownames(uv_summ), 2),
                         sep = "_")
  
  return(res)
})) %>% as.data.frame()

GGally::ggpairs(wind_agg[, grep("direction", names(wind_agg))])
GGally::ggpairs(wind_agg[, grep("speed", names(wind_agg))])
# Use mean, because it's the one with less eastern winds
# (although it's similar to median)

# Export data -------------------------------------------------------------

# merge with FWI data
fwi_data <- read.csv("data/patagonian_fires_data_with_fwi.csv")
clim_data <- cbind(fwi_data, wind_agg)

# add column to edit winddir and use in windninja
clim_data$direction_use <- wind_agg$direction_mean
clim_data$direction_editted <- 0
# write.csv(clim_data, "data/climatic_data_by_fire_FWI-wind.csv",
#           row.names = F)

# Some fires were separated for fire spread, and they are together in this database.
# caution when merging.