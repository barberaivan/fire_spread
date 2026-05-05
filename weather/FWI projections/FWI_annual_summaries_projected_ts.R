library(terra)
library(ggplot2)
library(viridis)
library(lubridate)
library(patchwork)
theme_set(theme_classic())

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),
    
    axis.line = element_line(linewidth = 0.35),
    
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# Temporal imputation of missing values across the time-series of a variable
impute_time <- function(x) {
  dseq <- 1:length(x)
  na_cols <- is.na(x)
  x[na_cols] <- approx(x = dseq[!na_cols], y = x[!na_cols],
                       xout = dseq[na_cols], rule = 2)$y
  return(x)
}

summarise <- function(x) {
  c(
    "n30" = sum(x > 30),
    "n38" = sum(x > 38),
    "p95" = quantile(x, 0.95, method = 8) |> unname()
  )
}

get_fseason <- function(d) {
  mm <- month(d)
  yy <- year(d)
  season <- yy
  season[mm >= 7] <- yy[mm >= 7] + 1
  return(season)
}

# Data --------------------------------------------------------------------

proj_dir <- file.path("data", "fwi_projections", "Quilcaille_Batibeniz_2023_database-nw_patagonia_clipped")
target_dir0 <- file.path("data", "fwi_projections")

# region of interest
# roi <- vect("/home/ivan/Insync/patagonian_fires/study_area/study_area.shp")
# roi <- vect(file.path(target_dir0, "points_for_fwi", "points_for_fwi.shp"))
roi <- vect(cbind(x = -71.5950379973438, y = -41.0017929442147), "points")
# Point in the PNNH, in forest. Using the polygon retrieved
# very smoothed FWI values, very low.

# Tables with climatic models and runs data
mtable <- read.csv(file.path(target_dir0, "models_table.csv"))
modmem <- read.csv(file.path(target_dir0, "models_emembers_table.csv"))
ff <- list.files(proj_dir)
scn <- c("historical", "ssp126", "ssp245", "ssp370", "ssp585")

modmem_ord <- modmem[order(modmem$model, modmem$emember), ]

# Select a single run for each model
mtable$emember <- "r1i1p1f1" # default
# As default is missing in some models, choose another.
sapply(mtable$model, function(x) {
  ("r1i1p1f1" %in% modmem$emember[modmem$model == x])
})
mtable$emember[mtable$model == "ACCESS-CM2"] <- "r2i1p1f1"
mtable$emember[mtable$model == "MIROC-ES2L"] <- "r1i1p1f2"

# Extract FWI values from one run by model --------------------------------

# By model and scenario, get the FWI ts in the roi
tslist <- vector("list", nrow(mtable))
names(tslist) <- mtable$model

for (i in 1:nrow(mtable)) {
  mod <- mtable$model[i]
  mem <- mtable$emember[i]
  
  message(paste(mod, "--------------------"))
  
  pattern <- paste0("(?=.*", mod, ")(?=.*", mem, ")")
  
  # List and filter files
  fflocal <- ff[grep(pattern, ff, perl = T)]
  
  # Loop over scenarios
  fwits <- lapply(1:5, function (j) {
    message(scn[j])
    # load raster
    rr <- rast(file.path(proj_dir, fflocal[j]))
    
    # Extract dates
    dates <- time(rr)
    
    # Extract values in roi
    fwi_raw <- terra::extract(rr, roi, method = "bilinear", raw = T)[, -1]
    
    # Impute missing values
    fwi_imputed <- impute_time(fwi_raw)
    
    # make df with all data.
    return(data.frame(date = dates, fwi = fwi_imputed))
  })
  
  tslist[[i]] <- fwits
}
saveRDS(tslist, file.path(target_dir0, "fwi_historical_proj_tslist.rds"))

# Summarise by year -------------------------------------------------------

tslist <- readRDS(file.path(target_dir0, "fwi_historical_proj_tslist.rds"))
nm <- nrow(mtable)

hbegin <- "1850-07-01"
hend <- "2014-06-30"
rows_hist <- date_historical >= hbegin & date_historical <= hend

summ_df <- do.call("rbind", lapply(1:nm, function(i) {
  print(paste("model", i))
  tss <- tslist[[i]]

  # Historical
  df_hist <- tss[[1]]
  df_hist$season <- get_fseason(df_hist$date)
  date_filter <- 
    df_hist$date >= hbegin & 
    df_hist$date <= hend
  hist_agg <- aggregate(fwi ~ season, df_hist[date_filter, ], summarise)
  hist_agg <- do.call("data.frame", hist_agg)
  hist_agg$scenario <- "historical"
  hist_agg$model <- mtable$model[i]
  hist_agg$emember <- mtable$emember[i]
  
  # Projections
  proj_agg <- do.call("rbind", lapply(2:5, function(j) {
    # print(j)
    # j = 2
    df_proj <- tss[[j]]
    df_proj$season <- get_fseason(df_proj$date)
    # merge historical ending with proj data
    dat <- rbind(
      df_hist[df_hist$date >= "2014-07-01", ],
      df_proj[df_proj$date <= "2100-06-30", ]
    )
    
    out_agg <- aggregate(fwi ~ season, dat, summarise)
    out_agg <- do.call("data.frame", out_agg)
    out_agg$scenario <- scn[j]
    out_agg$model <- mtable$model[i]
    out_agg$emember <- mtable$emember[i]
    
    return(out_agg)
  }))
  
  summ_list[[i]] <- rbind(hist_agg, proj_agg)
}))

# Average across models
summ_df_avg <- aggregate(
  cbind(fwi.n30, fwi.n38, fwi.p95) ~ season + scenario,
  summ_df, mean
)

saveRDS(summ_df, file.path(target_dir0, "fwi_summ_df.rds"))
saveRDS(summ_df_avg, file.path(target_dir0, "fwi_summ_df_avg.rds"))

# Plots -------------------------------------------------------------------

summ_df <- readRDS(file.path(target_dir0, "fwi_summ_df.rds"))
summ_df_avg <- readRDS(file.path(target_dir0, "fwi_summ_df_avg.rds"))

summ_df <- summ_df[summ_df$season >= 1950, ]
summ_df_avg <- summ_df_avg[summ_df_avg$season >= 1950, ]

# Average
a <-
ggplot(summ_df_avg, aes(season, fwi.n30, color = scenario)) + 
  geom_line() + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("Days with FWI > 30") + 
  nice_theme() +
  theme(legend.title = element_blank())
  
b <-
ggplot(summ_df_avg, aes(season, fwi.n38, color = scenario)) + 
  geom_line() + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("Days with FWI > 38") +
  theme(legend.title = element_blank()) + 
  nice_theme()

c <-
ggplot(summ_df_avg, aes(season, fwi.p95, color = scenario)) + 
  geom_line() + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("FWI percentile 95 %") +
  theme(legend.title = element_blank()) +
  nice_theme()
  

merged1 <- (
  a + theme(legend.position = "none",
            axis.title.x = element_blank()) +
  b + theme(legend.position = "bottom") +
  c + theme(legend.position = "none",
            axis.title.x = element_blank())
) 
ggsave(file.path(target_dir0, "FWI_projections_ts.png"), plot = merged1,
       width = 16, height = 9, units = "cm", dpi = 900)

ggsave(file.path(target_dir0, "FWI_projections_ts_gt30.png"), plot = a,
       width = 11, height = 9, units = "cm", dpi = 900)
ggsave(file.path(target_dir0, "FWI_projections_ts_gt38.png"), plot = b,
       width = 11, height = 9, units = "cm", dpi = 900)
ggsave(file.path(target_dir0, "FWI_projections_ts_p95.png"), plot = c,
       width = 11, height = 9, units = "cm", dpi = 900)

# Average + realizations
A <- 
ggplot() + 
  geom_line(data = summ_df,
            mapping = aes(season, fwi.n30, color = scenario, group = model),
            alpha = 0.3, linewidth = 0.4) + 
  geom_line(data = summ_df_avg,
            mapping = aes(season, fwi.n30, color = scenario),
            linewidth = 0.8) + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("Days with FWI > 30")

B <- 
ggplot() + 
  geom_line(data = summ_df,
            mapping = aes(season, fwi.n38, color = scenario, group = model),
            alpha = 0.3, linewidth = 0.4) + 
  geom_line(data = summ_df_avg,
            mapping = aes(season, fwi.n38, color = scenario),
            linewidth = 0.8) +  
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("Days with FWI > 38")  

C <- 
ggplot() + 
  geom_line(data = summ_df,
            mapping = aes(season, fwi.p95, color = scenario, group = model),
            alpha = 0.3, linewidth = 0.4) + 
  geom_line(data = summ_df_avg,
            mapping = aes(season, fwi.p95, color = scenario),
            linewidth = 0.8) + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Year") +
  ylab("FWI percentile 95 %")  


A+B+C

# Export table ------------------------------------------------------------

write.csv(summ_df_avg, file.path(target_dir0, "fwi_annual_summaries.csv"),
          row.names = F)



# Perc95 in Spanish, for thesis defense -----------------------------------

d2 <- summ_df_avg
d2$scenario <- factor(
  d2$scenario, 
  levels = unique(d2$scenario),
  labels = c("Histórico", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
)

ggplot(d2, aes(season, fwi.p95, color = scenario)) + 
  geom_line() + 
  scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  xlab("Año") +
  ylab("FWI (percentil 95 %)") +
  theme(legend.title = element_blank()) +
  nice_theme()


# Get modern record -------------------------------------------------------

# Import image
rr <- rast("data/fwi_daily_1998-2022/24km/fwi_daily_19970701_20230630.tif")

# Extract dates
dates <- time(rr)

# Extract values in roi
fwi_raw <- terra::extract(rr, roi, method = "bilinear", raw = T)[, -1]

# Impute missing values
fwi_imputed <- impute_time(fwi_raw)

# make df with all data.
dmod <- data.frame(date = dates, fwi = fwi_imputed)

# Get fire season
dmod$season <- get_fseason(dmod$date)

# Aggregate (perc95 by season)
dmod <- aggregate(fwi ~ season, dmod, quantile, prob = 0.95, method = 8)

# Merge with historical scale
dhist <- d2[d2$scenario == "Histórico", ]
dhist_comp <- dhist[dhist$season %in% out_agg$season, ]
dmod_comp <- dmod[dmod$season %in% dhist$season, ]

# scale to historical range
dmod2 <- dmod
fwiz <- (dmod$fwi - mean(dmod_comp$fwi)) / sd(dmod_comp$fwi)
dmod2$fwi <- fwiz * sd(dhist_comp$fwi.p95) + mean(dhist_comp$fwi.p95)

# Merge with d2
names(d2)
d3 <- d2[, c("season", "scenario", "fwi.p95")]
names(d3) <- c("season", "scenario", "fwi")
d3$scenario <- as.character(d3$scenario)

dmod2$scenario <- "Moderno"
dmod2 <- dmod2[, c("season", "scenario", "fwi")]

d4 <- rbind(d3, dmod2)
d4$scenario <- factor(
  as.character(d4$scenario), 
  levels = c("Histórico", "Moderno", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
)

# Ribbons for periods
ylow <- 25; yhigh <- 44
ribs <- data.frame(
  season = c(1999, 2022, 2040, 2049, 2090, 2099),
  period = rep(c("modern", "40s", "90s"), each = 2),
  ylow = rep(ylow, 6),
  yhigh = rep(yhigh, 6)
)


filter <- 
  d4$scenario != "Histórico" &
  (d4$scenario != "Moderno" & d4$season >= 2023) |
  (d4$scenario == "Moderno" & d4$season <= 2023) 

ggplot(
  data = d4[filter, ], 
  mapping = aes(season, fwi, color = scenario)
  ) + 
  geom_ribbon(
    data = ribs, mapping = aes(season, ymin = ylow, ymax = yhigh, group = period),
    inherit.aes = F, alpha = 0.15, color = NA
  ) +
  geom_line(alpha = 0.8) + 
  geom_smooth(se = F, method = "loess", span = 1) +
  # scale_color_viridis(discrete = T, end = 0.9, option = "A") +
  scale_color_viridis(option = "C", discrete = T, end = 0.8, 
                      name = "Escenario\nclimático") +
  xlab("Año") +
  ylab("FWI (percentil 95 %)") +
  theme(legend.title = element_blank()) +
  nice_theme() + 
  scale_y_continuous(limits = c(ylow, yhigh), expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(
  "fire regime simulations/figures_defensa/FWI_projections.png",
  width = 14, height = 10, units = "cm", dpi = 600
)