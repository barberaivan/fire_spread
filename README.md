# Fire spread and ignitions modelling

Fitting fire spread and ignition models from fire maps, using Approximate Bayesian Computation.

The fire spread model is a cellular automata where fire can spread in the 8-neighbourhood of each burning cell ("queen" direction following the terra R package nomenclature). The spread probability follows a Bernoully distribution, which probability is defined as a linear model at the logit scale. The spread probability varies a function of landscape variables that can vary in space and time (climatic interannual variation).
The rasters are treated as the terra R package does: vectors with a cell id for every pixel. The landscape where fire can spread is a raster stack with a layer by covariate (approximately), but in the vectorial representation it is stored as a matrix with a row by pixel and a column by layer.
As simulating spread is expensive we write the simulation function in C++ and use the Rcpp R package to make the function available in R.  
  
We use large files storing the landscape data, so we don't upload them to the repo. The data can be asked to ivanbarbera93@gmail.com, and the fire_spread_data folder must be placed in the same directory as the cloned repo to work (e.g., ```/a_path/fire_spread``` and ```/a_path/fire_spread_data```). The working directory has to be set at ```/a_path/fire_spread```.