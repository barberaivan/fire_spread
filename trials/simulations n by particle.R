fire_sim <- readRDS("data/simulations/fire_sim.rds")
str(fire_sim)
length(fire_sim) # 200 fuegos simulados

# array con las metricas de comparacion
fire_sim_disc_arr <- readRDS("data/simulations/fire_sim_disc_arr.rds")

str(fire_sim_disc_arr)
ovmat <- fire_sim_disc_arr[, , "overlap_sp"]
View(ovmat)
