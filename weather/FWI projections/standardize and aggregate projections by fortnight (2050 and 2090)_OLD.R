# Packages ----------------------------------------------------------------

library(terra)
source(file.path("weather", "fortnight_functions.R"))


# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# Bhattacharyya distance between empirical densities, computed from samples x and y.
# n is the number of points to compute the initial density, then it is
# multiplied by 3. lower is the lower limit to fit the density.
# wx and wy are the observation weights for x and y.
bhattacharyya_distance <- function(x, y, wx, wy, n = 2 ^ 11, lower = NULL,
                                   similarity = TRUE) {
  # Initial densities, to check limits
  if(is.null(lower)) {
    dx0 <- density(x, weights = wx, n = n)
    dy0 <- density(y, weights = wy, n = n)
  } else {
    dx0 <- density(x, weights = wx, from = lower, n = n)
    dy0 <- density(y, weights = wy, from = lower, n = n)
  }
  
  if(is.null(lower)) {
    lower <- min(c(dx0$x, dy0$x))
  }
  upper <- max(c(dx0$x, dy0$x))
  # refit densities with lower and upper
  dx <- density(x, weights = wx, from = lower, to = upper, n = n * 3)
  dy <- density(y, weights = wy, from = lower, to = upper, n = n * 3)
  
  sqrt_prod <- function(x) {
    dens_x <- approx(dx$x, dx$y, xout = x, rule = 2)$y
    dens_y <- approx(dy$x, dy$y, xout = x, rule = 2)$y
    return(sqrt(dens_x * dens_y))
  }
  
  bc <- integrate(sqrt_prod, lower, upper, subdivisions = 1e5)$value
  
  if(similarity) return(bc) else 
  return(-log(bc))
}

# Temporal imputation of missing values across the cells of a raster.
raster_impute_time <- function(r) {
  va <- values(r)
  dseq <- 1:ncol(va)
  for(p in 1:ncell(r)) {
    na_cols <- is.na(va[p, ])
    va[p, na_cols] <- approx(x = dseq[!na_cols], y = va[p, !na_cols],
                             xout = dseq[na_cols], rule = 2)$y
  }
  values(r) <- va
  return(r)
}

# Aggregate raster by fortnights. It must have a time property to compute
# fortnights.
aggregate_fortnights <- function(r) {
  # extract date and fortnight
  dates <- time(r)
  forts_all <- date2fort(dates)
  forts <- unique(forts_all)
  
  # loop to average
  r_fort <- do.call("c", lapply(forts, function(f) {
    ids <- which(forts_all == f)
    return(mean(r[[ids]]))
  })) 
  
  # set new time and name.
  time(r_fort) <- fort2date(forts)
  names(r_fort) <- forts
  
  return(r_fort)
}


# Data --------------------------------------------------------------------

# pnnh vector to crop rasters
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]

# directories for present data (ERA5, 24km) and projections
data_dir <- file.path("data", "fwi_daily_1998-2022", "24km")
proj_dir <- file.path("data", "fwi_projections", "Quilcaille_Batibeniz_2023_database-nw_patagonia_clipped")
target_dir <- file.path("data", "fwi_projections", "fwi_fortnights_standardized")
target_dir2 <- file.path("data", "fwi_projections", "fwi_fortnights_standardized_modern_compare")
target_dir0 <- file.path("data", "fwi_projections")

proj_files <- list.files(proj_dir)

# modern data raster
modern_daily <- rast(file.path(data_dir, "fwi_daily_19970701_20230630.tif"))
modern_fort <- rast(file.path(data_dir, "fwi_fortnights_19970701_20230630_standardized.tif"))

# To compare values in the PNNH, project to posgar 2007-1 and crop to the park.
r0 <- project(modern_fort, "EPSG:5343", method = "cubicspline", threads = 8)
pnnhf <- project(pnnh, "EPSG:5343")
raster_ref <- crop(r0, pnnhf, snap = "out")

# write this cropped and projected raster to file
# writeRaster(raster_ref, file.path(data_dir, "fwi_fortnights_19970701_20230630_standardized_pnnh.tif"))

# plot(raster_ref)
# plot(pnnhf, add = T)
# res(raster_ref)
# ncell(raster_ref)
vobs0 <- terra::extract(raster_ref, pnnhf, exact = T, raw = T) # values with pixel weights
out <- which(colnames(vobs0) %in% c("ID", "fraction"))
vobs_mat <- vobs0[, -out]
wobs <- normalize(vobs0[, "fraction"])

# Models-ensemble members available at all scenarios -----------------------

proj_files
nf <- length(proj_files)

projdata <- data.frame(
  model = character(nf),
  scenario = character(nf),
  emember = character(nf),
  filename = proj_files
)

for(i in 1:nf) {
  projdata[i, 1:3] <- strsplit(proj_files[i], 
                               split = "fwi_day_|_|_native.nc")[[1]][2:4]
}

models <- unique(projdata$model)
nm <- length(models)

scenarios <- unique(projdata$scenario)
ns <- length(scenarios)

emembers <- unique(projdata$emember)
ne <- length(emembers)

# combinations of models-ensemble members available
modmem_data <- projdata[!duplicated(projdata[, c("model", "emember")]), ]

# models-ensemble members available for all scenarios, including historical
modmem_data$complete <- F
for(i in 1:nrow(modmem_data)) {
  filt <- 
    projdata$model == modmem_data$model[i] &
    projdata$emember == modmem_data$emember[i]
  scn_available <- projdata$scenario[filt]
  modmem_data$complete[i] <- all(scenarios %in% scn_available)
}
sum(modmem_data$complete) # 157 (*4).

models_n <- aggregate(emember ~ model, modmem_data[modmem_data$complete, ], 
                      length)
models_n2 <- aggregate(emember ~ model, modmem_data, length)

# which have no crs?

models_crs <- models_n2
models_crs$crs <- NA 

emat <- matrix(NA, nrow(models_crs), 4)
colnames(emat) <- c("xmin", "xmax", "ymin", "ymax")

resmat <- matrix(NA, nrow(models_crs), 2)
colnames(resmat) <- c("x", "y")

rcmat <- matrix(NA, nrow(models_crs), 2)
colnames(rcmat) <- c("rows", "cols")

# for(m in 1:nrow(models_crs)) {
#   # m = 1
#   fn <- proj_files[grep(models_crs$model[m], proj_files)][1]
#   rr <- rast(file.path(proj_dir, fn))
#   
#   emat[m, ] <- as.vector(ext(rr))
#   models_crs$crs[m] <- crs(rr)
#   resmat[m, ] <- res(rr)
#   rcmat[m, 1] <- nrow(rr)
#   rcmat[m, 2] <- ncol(rr)
# }

models_crs$ext <- emat
models_crs$res <- resmat
models_crs$dim <- rcmat

write.csv(models_crs, file.path(target_dir0, "models_crs.csv"),
          row.names = F)



# Compute fwi fort z for all periods, models, members and scenarios -------

modmem <- modmem_data[modmem_data$complete, c("model", "emember")]
nn <- nrow(modmem)
rownames(modmem) <- 1:nn

nsp <- ns - 1 # number of projection scenarios

# begin and end for each period. The reference period is the one used for 
# standardization, but also to compare modern climate reanalysis with 
# simulations. Note that this period takes the historical and the beggining 
# of climatic scenarios.
bref <- as.Date("1998-07-01", format = "%Y-%m-%d")
eref <- as.Date("2023-06-30", format = "%Y-%m-%d")

hist_end <- as.Date("2014-12-31", format = "%Y-%m-%d")
proj_begin <- as.Date("2015-01-01", format = "%Y-%m-%d")

# begin and end of target decades
bdec <- as.Date(c("2039-07-01", "2089-07-01"), format = "%Y-%m-%d")
edec <- as.Date(c("2049-06-30", "2099-06-30"), format = "%Y-%m-%d")

dec_names <- c("2040", "2090")
  
# tables to save na proportion in the daily database in projection ts
naprop <- list(
  matrix(NA, nrow(modmem), nsp),
  matrix(NA, nrow(modmem), nsp)
)
colnames(naprop[[1]]) <- colnames(naprop[[2]]) <- scenarios[-1]

# table to save the bhattacharyya similarity between projected and reanalized
# reference period.
batsim <- naprop[[1]]

# subset dates to compare in modern climate
colcol <- as.numeric(colnames(vobs_mat))
cols_comp <- colcol >= date2fort(bref) & colcol <= date2fort(eref)
vobs_mat_comp <- vobs_mat[, cols_comp]

vobs_comp <- as.vector(vobs_mat_comp)  # all pixel values, but weighted by
wobs_comp <- normalize(rep(wobs, sum(cols_comp))) # representativity in PNNH

# modmem$crs <- NA
# modmem$ext <- NA
# extmat <- naprop40
# colnames(extmat) <- c("xmin", "xmax", "ymin", "ymax")
# 
# for(i in 1:nn){
#   # i = 1
# mod <- modmem$model[i]
# mem <- modmem$emember[i]
# 
# # Load rasters for the 5 scenarios
# filt <- projdata$model == mod & projdata$emember == mem
# filenames <- projdata$filename[filt]
# 
# # historical
# shist <- rast(file.path(proj_dir, filenames[1]))
# # crs(shist)
# # library(ncdf4)
# # zzz <- nc_open(file.path(proj_dir, filenames[1]))
# # crds(shist)
# modmem$crs[i] <- crs(shist)
# ee <- ext(shist)
# extmat[i, ] <- as.vector(ee)
# }
# modmem$ext <- extmat
# write.csv(modmem, file.path(target_dir2, "tabla.csv"))
tabla <- read.csv(file.path(target_dir0, "tabla.csv"))
i_ok <- which(!tabla$ext.xmin == 0)

# t1 <- tabla[tabla$ext.xmin != 0, ]
# t2 <- tabla[!duplicated(tabla[, "model"]), ]
# t2$has_CRS <- 0
# t2$has_CRS[t2$ext.xmin != 0] <- 1
# write.csv(t2[, c("model", "has_CRS")], file.path(target_dir0, "tabla2.csv"),
#           row.names = F)
# 

t3 <- tabla[!duplicated(tabla[, "model"]), ]
t3$has_CRS <- 0
t3$has_CRS[t3$ext.xmin != 0] <- 1
t3 <- t3[, grep("model|has_CRS|ext.", colnames(t3))]
write.csv(t3, file.path(target_dir0, "tabla2_models.csv"),
          row.names = F)

## IMPORTANTE
# muchos modelos no tienen CRS. Hasta corregir eso, uso los que sí tienen.
# la primera corrida hizo los primeros 5 sin drama, porque tenían.
# corro los i_ok > 5 por ahora

# Loop over model-members
for(i in i_ok[i_ok>=96]) { ## for(i in 2:nn) {
  # i = i_ok[1]
  mm <- paste(i, ":", modmem$model[i], modmem$emember[i],
              "---------------------------------")
  message(mm)

  mod <- modmem$model[i]
  mem <- modmem$emember[i]
  
  filt <- projdata$model == mod & projdata$emember == mem
  filenames <- projdata$filename[filt]
  
  # historical
  rhist <- rast(file.path(proj_dir, filenames[1]))
  layers_hist <- which(time(rhist) >= bref)
  rhist <- rhist[[layers_hist]]
  
  # loop over each projection scenario
  for(s in 1:nsp) {
    # s = 1
    message(scenarios[s+1])
    
    rproj0 <- rast(file.path(proj_dir, filenames[s+1]))
    
    # complete the reference period with the projection
    layers_ref_end <- which(time(rproj0) <= eref)
    rref <- c(rhist, rproj0[[layers_ref_end]])
    
    # fill NA in reference period with linear interpolation in time
    rref <- raster_impute_time(rref)
    
    # constants to standardize
    rref_mean <- mean(rref, na.rm = T)
    rref_sd <- app(rref, sd, na.rm = T)
    
    # standardize reference period (for comparison with modern climate)
    rrefz <- (rref - rref_mean) / rref_sd
    
    # aggregate by fortnight
    rrefz_fort <- aggregate_fortnights(rrefz)
    
    # reproject to posgar in PNNH (final step)
    rref_out <- project(rrefz_fort, raster_ref, method = "cubicspline", 
                        threads = 8)
    
    # write to disk
    nncomp <- paste("fwi-fortnight-z",  
                    mod, scenarios[s+1], mem, 
                    "modern-compare", "pnnh.tif", sep = "_")
    writeRaster(rref_out, file.path(target_dir2, nncomp), overwrite = T)
    
    # compare distribution 
    vhist0 <- terra::extract(rref_out, pnnhf, exact = T, raw = T)
    out <- which(colnames(vhist0) %in% c("ID", "fraction"))
    vhist_mat <- vhist0[, -out]
    whist <- normalize(rep(vhist0[, "fraction"], ncol(vhist_mat)))
    vhist <- as.vector(vhist_mat)
    
    # Compute bat
    batsim[i, s] <- bhattacharyya_distance(vhist, vobs_comp,
                                           whist, wobs_comp)
    

    # Prepare raster for the two target decades __________________________
    for(d in 1:2) {
      
      # filter projection periods and standardize
      layers_focal <- which(time(rproj0) >= bdec[d] & time(rproj0) <= edec[d])
      rproj <- rproj0[[layers_focal]]
      rprojz <- (rproj - rref_mean) / rref_sd
      
      # get amount of daily NA values
      vvv <- as.vector(values(rprojz))
      naprop[[d]][i, s] <- sum(is.na(vvv)) / length(vvv)
      
      # fill NA with linear interpolation in time
      rprojz <- raster_impute_time(rprojz)
      
      # aggregate by fortnight
      rprojz_fort <- aggregate_fortnights(rprojz)
      
      # resample to PNNH extent, and reproject to posgar 2007-1
      rproj_out <- project(rprojz_fort, raster_ref, method = "cubicspline", 
                           threads = 8)
      
      # write results
      nnn <- paste("fwi-fortnight-z", mod, scenarios[s+1], mem,  
                   dec_names[d], "pnnh.tif", sep = "_")
      writeRaster(rproj_out, file.path(target_dir, nnn), overwrite = T)
    } # end loop over decades
  } # end loop over scenarios
  
  # save tables if 10 iters have passed
  if(i %% 10 == 0) {
    write.csv(batsim, file.path("data", "fwi_projections", "temptable_batsim.csv"))
    write.csv(naprop[[1]], file.path("data", "fwi_projections", "temptable_naprop40.csv"))
    write.csv(naprop[[2]], file.path("data", "fwi_projections", "temptable_naprop90.csv"))
  }
}

write.csv(batsim, file.path("data", "fwi_projections", "temptable_batsim.csv"))
write.csv(naprop[[1]], file.path("data", "fwi_projections", "temptable_naprop40.csv"))
write.csv(naprop[[2]], file.path("data", "fwi_projections", "temptable_naprop90.csv"))
# done




# CLEAN ALL WHEN WE HAVE ALL MODELS WITH CRS ------------------------------



