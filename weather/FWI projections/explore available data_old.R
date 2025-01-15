# Packages ----------------------------------------------------------------

library(terra)

# Data --------------------------------------------------------------------

# directories for present data (ERA5, 24km) and projections
data_dir <- file.path("data", "fwi_daily_1998-2022", "24km")
proj_dir <- file.path("data", "fwi_projections", "Quilcaille_Batibeniz_2023_database-nw_patagonia_clipped")
proj_files <- list.files(proj_dir)

# modern data raster
fwi_modern_full <- rast(file.path(data_dir, "fwi_daily_19970701_20230630.tif"))

# time range for comparison/calibration
begin <- as.Date("1998-01-01")
end <- as.Date("2014-12-31")

timefilt1 <- 
  time(fwi_modern_full) >= begin &
  time(fwi_modern_full) <= end
fwi_modern <- fwi_modern_full[[which(timefilt1)]]

# Example raster for historical scenario
filename <- file.path(proj_dir, "fwi_day_MIROC6_historical_r10i1p1f1_native.nc")
r1 <- rast(filename)
timefilt <- time(r1) >= begin & time(r1) <= end
r1 <- r1[[which(timefilt)]]


r1_resample <- project(r1, fwi_modern)

vobs <- as.vector(values(fwi_modern))
vhist <- as.vector(values(r1_resample))
mm <- lm(vhist ~ vobs)
summary(mm)

plot(fwi_modern[[1]])
plot(r1_resample[[1]])

hist(vobs[!is.na(vhist)])
hist(vhist)

# resample from observed to historical
obs_resample <- project(fwi_modern, r1)

plot(mean(obs_resample))
plot(mean(r1, na.rm = T))

vobs <- values(obs_resample)
vhist <- values(r1)

par(mfrow = c(2, 1))
p = 3
plot(vobs[p, ], type = "l")
plot(vhist[p, ], type = "l")
par(mfrow = c(1, 1))
# Están en escalas muy distintas

plot(vobs[p, ] ~ vhist[p, ], col = rgb(0 ,0, 0, 0.05))
plot(vhist[p, ], type = "l")


# Hay discrepancias muy grandes, no solo en la escala sino también en el 
# promedio a lo largo del tiempo. 
# Hacer una lista de todos los modelos que hay y su resolución,
# a ver si encontramos algo a 0.5 °

# Otra opción es cuantificar una distribución de aumentos a nivel anual 
# (tipo mean futuro / mean historical)
# y usar la distribución presente multiplicada por esa distribución
# de factores.

# O agarrar cada modelo-ensamble, estandarizar en base al período 1998-2022
# Y luego usar esas anomalías



# Proyecciones disponibles ------------------------------------------------

proj_files
nf <- length(proj_files)

projdata <- data.frame(
  model = character(nf),
  scenario = character(nf),
  ensemble = character(nf)
)

for(i in 1:nf) {
  projdata[i, ] <- strsplit(proj_files[i], 
                            split = "fwi_day_|_|_native.nc")[[1]][2:4]
}
View(projdata)

models <- unique(projdata$model)
nm <- length(models)

scenarios <- unique(projdata$scenario)
ns <- length(scenarios)

ensembles <- unique(projdata$ensemble)
ne <- length(ensembles)

projdata2 <- aggregate(ensemble ~ model + scenario, projdata, length)
projdata3 <- aggregate(ensemble ~ model, projdata, length)
projdata4 <- projdata[!duplicated(projdata[, c("model", "ensemble")]), ]

# Which scenarios are available for each model?
scn <- matrix(0, nrow(projdata3), ns)
for(m in 1:nm) {
  for(s in 1:ns) {
    filt <- 
      projdata$model == projdata3$model[m] & 
      projdata$scenario == scenarios[s]
    scn[m, s] <- nrow(projdata[filt, ])
  }
}

colnames(scn) <- scenarios
scn <- as.data.frame(scn)

modeldata <- cbind(projdata3[, 1, drop = F], scn)

# get spatial resolution for each model

modeldata$res_x <- modeldata$res_y <- numeric(nm)

for(m in 1:nm) {
  print(m)
  mod <- modeldata$model[m]
  fileget <- proj_files[grep(mod, proj_files)[1]]
  rr <- rast(file.path(proj_dir, fileget))
  resol <- res(rr)
  modeldata$res_x[m] <- resol[1]
  modeldata$res_y[m] <- resol[2]
}

write.csv(
  modeldata, 
  file.path("data", "fwi_projections", "fwi_projections_models_ensemble-count.csv")
)

# Do all ensembles for each scenario have the equivalent ensemble in the 
# historical?

ematch <- modeldata[, c("model", scenarios[2:ns])]
ematch[, 2:5] <- NA
for(m in 1:nm) {
  for(s in 2:ns) {
    
    if(modeldata[m, s+1] == 0) next
      
    filt1 <-
      projdata$model == modeldata$model[m] &
      projdata$scenario == scenarios[s]
    ensembles_sc <- unique(projdata$ensemble[filt1])
    
    filt2 <-
      projdata$model == modeldata$model[m] &
      projdata$scenario == "historical"
    ensembles_hist <- unique(projdata$ensemble[filt2])
    
    ematch[m, s] <- sum(ensembles_sc %in% ensembles_hist)
  }
}

ematch
ematch2 <- cbind(ematch, modeldata[, c("res_x", "res_y")])
write.csv(
  ematch2, 
  file.path("data", "fwi_projections", "fwi_projections_models_ensemble-count-match-historical.csv")
)



# Probamos EC-Earth3

# Se puede unir la serie de historical con ssp1?
fproj <- proj_files[grep("fwi_day_EC-Earth3_ssp126_", proj_files)][1]
fhist <- proj_files[grep("fwi_day_EC-Earth3_historical_", proj_files)][2]

# same ensemble, contrasting scn
rh <- rast(file.path(proj_dir, fhist))
rp <- rast(file.path(proj_dir, fproj))
range(time(rh))
range(time(rp))

# Pasos:

# Usar sólo modelos-ensambles presentes en los 4 escenarios, para que en todos
# se hayan usado los mismos.

# models-ensembles available:
projdata4
# identify those present in all scenarios
projdata4$complete <- F
for(i in 1:nrow(projdata4)) {
  filt <- 
    projdata$model == projdata4$model[i] &
    projdata$ensemble == projdata4$ensemble[i]
  scn_available <- projdata$scenario[filt]
  projdata4$complete[i] <- all(scenarios %in% scn_available)
}
sum(projdata4$complete) # 157!! a lot

# Parear histórico-scenario por modelo y ensamble,
# filtrar time periods de interés [97-22, 51-60, 91-100],
# estandarizar igual que con clima moderno, 
# promediar por quincena, 

# resamplear a escala de interés.
# guardar




# Ejemplo con 1 modelo-ensamblemember -------------------------------------










