# Functions and objects to use fortnights as temporal unit.

library(lubridate)
library(dplyr)

# Constants/Data ----------------------------------------------------------

# Create a reference table with date and fortnight-date, to match observations.
# The fortnight-date (14 days) is labelled as ending-year_fortnight.
# Years run from July to June, named by the year of the ending month. Fortnights
# are centred at January first, so irregular-length ones occur at the first or
# last month, which correspond to the southern winter (less relevant for fire).
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
# plot(dtable$fort)
# max(dtable$fort) / length(years) # OK

# THE CEMS FWI WAS COMPUTED WITH THESE FUNCTIONS, SO ORIGIN MUST BE 1996
# (FORTNIGHT 1 BEGINS IN 1995-07-01).

# Functions ---------------------------------------------------------------

# based on the table, turn date into (continuous) fortnight
date2fort <- function(d) {
  dd <- data.frame(date = d)
  dd <- left_join(dd, dtable, by = "date")
  if(anyNA(dd$fort)) {
    warning("Found dates outside the reference table, returning NA.")
  }
  return(dd$fort)
}

# based on the table, turn date into (continuous) fortnight

# number of fortnights of the upper half of the year (january-june)
days_below <- as.numeric(as.Date("1998-12-31") - as.Date("1998-07-01") + 1)
forts_half <- round(days_below / 14, 0) |> as.integer() # number of fortnights by half year
nfy <- forts_half * 2 # number of fortnights by year

 
date2fort <- function(d, origin = 1997) {
  # d <- as.Date(c("1998-07-01", "1999-02-13", "2000-04-30", "2099-06-30", "2000-08-04"))
  # d <- as.Date(c("1998-07-01"))
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
    if(is.Date(origin)) origin <- year(origin)
    out <- out + nfy * (y - origin)
  }
  
  return(out)
}

yy <- 2100
date2fort(paste(yy, "01-06", sep = "-"), origin = yy)


# based on the table, turn fortnight (continuous) into the first date
fort2date <- function(f) {
  # f = c(-60, 120, 133) # test
  f <- sort(f)
  dates <- as.Date(character(length(f)), format = "%Y-%m-%d")
  use <- f %in% dtable$fort
  ids_use <- which(use)

  dat <- dtable[dtable$fort %in% f[use], ]
  dat <- dat[!duplicated(dat[, "fort"]), ] # dtable is ordered by date,
  # so it takes the first
  dates[use] <- dat$date
  if(anyNA(dates)) {
    warning("Found fortnights outside the reference table, returning NA.")
  }
  return(dates)
}