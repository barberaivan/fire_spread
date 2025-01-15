# Functions and objects to use fortnights as temporal unit.

# Create a reference table with date and fortnight-date, to match observations.
# The fortnight-date (14 days) is labelled as ending-year_fortnight.
# Years run from July to June, named by the year of the ending month. Fortnights
# are centred at January first, so irregular-length ones occur at the first or
# last month, which correspond to the southern winter (less relevant for fire).

# IMPORTANT: the fwi for the modern period was averaged by fortnight using 
# functions with non-flexible origin. So, for compatibility, the origin year
# must be set to 1996 (fortnight 1 starts at 1995-07-01).


# Packages ----------------------------------------------------------------

library(lubridate)
library(dplyr)

# Constants ---------------------------------------------------------------

days_below <- as.numeric(as.Date("1998-12-31") - as.Date("1998-07-01") + 1)
forts_half <- round(days_below / 14, 0) |> as.integer() # number of fortnights by half year
nfy <- forts_half * 2 # number of fortnights by year
origin_default <- 1996

# Functions ---------------------------------------------------------------

# origin must be a fire-year, running from July to June, labelled with the 
# year of the ending date
date2fort <- function(d, origin = origin_default) {
  
  d <- as.Date(d)
  m <- month(d)
  y <- year(d)
  
  out <- ifelse(
    m < 7,
    
    # ending half of the year (after january)
    ceiling( 
      as.numeric(d - as.Date(paste(y, "01-01", sep = "-")) + 1) / 14
    ) + forts_half, 
    
    # beggining half of the year
    forts_half - 
      pmin(
        ceiling(as.numeric(as.Date(paste(y, "12-31", sep = "-")) - d + 1) / 14),
        forts_half 
      ) + 1
  )
  
  if(!is.null(origin)) {
    fire_year <- ifelse(m > 6, y + 1, y)
    out <- out + nfy * (fire_year - origin)
  }
  
  return(out)
}

# get the first date of the fortnight
fort2date0 <- function(f, origin = origin_default) {
  # get fire_year
  year_rank <- ceiling(f / nfy)
  y <- origin + year_rank - 1
  
  # make a date sequence for the whole year and match
  # (brute, but easy to program)
  dfrom <- as.Date(paste(y-1, "07-01", sep = "-"))
  dto <- as.Date(paste(y, "06-30", sep = "-"))
  dateseq <- seq(dfrom, dto, 1)
  forts <- date2fort(dateseq, origin)
  date <- min(dateseq[forts == f])
  
  return(date)
}

# the same but looping
fort2date <- function(f, origin = origin_default) {
  dates <- Date(length(f))
  for(i in 1:length(f)) {
    dates[i] <- fort2date0(f[i], origin)
  }
  return(dates)
}

# Get all the fortnights for a given year
forts_from_year <- function(year, origin = origin_default) {
  start <- as.Date(paste(as.character(year-1), "07", "01", sep = "-"))
  end <- as.Date(paste(as.character(year), "06", "30", sep = "-"))
  seqd <- seq(start, end, 10)
  forts <- unique(date2fort(seqd))
  return(forts)
}

# Tests -------------------------------------------------------------------

# yy <- 2099
# dd <- as.Date(paste(yy, "06-25", sep = "-"))
# (ff <- date2fort(dd, origin = 2099))
# fort2date(ff, origin = 2099)
# 


# Objects for backwards compatibility -------------------------------------

nfy <- floor(365/14)
year_low <- 1996; year_high <- 2023
years <- year_low:year_high

# Date table
dtable <- do.call("rbind", lapply(years, function(y) {
  # y = 1997
  textlow <- paste(y-1, "07", "01", sep = "-")
  textmid <- paste(y, "01", "01", sep = "-")
  texthigh <- paste(y, "06", "30", sep = "-")
  
  dlow <- as.Date(textlow, format = "%Y-%m-%d")
  dmid <- as.Date(textmid, format = "%Y-%m-%d")
  dhigh <- as.Date(texthigh, format = "%Y-%m-%d")
  
  # upper half of the year
  nd_upper <- as.numeric(dhigh - dmid) + 1
  nf_upper <- floor(nd_upper / 14)
  rem_upper <- nd_upper %% 14
  
  if(rem_upper >= 7) {
    fort_upper <- rep(1:(nf_upper + 1), c(rep(14, nf_upper), rem_upper))
  } else {
    lengths <- rep(14, nf_upper)
    lengths[nf_upper] <- lengths[nf_upper] + rem_upper
    fort_upper <- rep(1:nf_upper, lengths)
  }
  
  df_up <- data.frame(date = seq(dmid, dhigh, 1),
                      fort0 = fort_upper)
  
  # lower half of the year
  nd_lower <- as.numeric(dmid - dlow)
  nf_lower <- floor(nd_lower / 14)
  rem_lower <- nd_lower %% 14
  
  if(rem_lower >= 7) {
    fort_lower <- rep(-nf_lower:0, c(rem_lower, rep(14, nf_lower)))
  } else {
    lengths <- rep(14, nf_lower)
    lengths[1] <- lengths[1] + rem_lower
    fort_lower <- rep(-(nf_lower-1):0, lengths)
  }
  
  df_low <- data.frame(date = seq(dlow, dmid - 1, 1),
                       fort0 = fort_lower)
  
  # merge and tidy
  dd <- rbind(df_low, df_up)
  dd$fort_focal <- dd$fort0 - min(dd$fort0) + 1
  dd$year <- as.numeric(y)
  
  return(dd[, c("date", "fort_focal", "year")])
}))

# Get continuous fortnight identifier, from year_low to year_high.
year_ref <- dtable$year - min(dtable$year)
dtable$fort <- dtable$fort_focal + year_ref * nfy
