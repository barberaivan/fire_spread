# Exploration of individual likelihoods.

# Here I explore the likelihoods approximated with GAMs, created in
# <fire_spread/pseudolikelihood estimation/pseudolikelihoods_estimation_ig-known.R>

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(ggdensity) # geom_hdr

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(posterior)     # manage posterior samples
library(tidybayes)     # not sure if this was used
library(bayesplot)     # visualize posteriors


# Functions ---------------------------------------------------------------

# the same as tidy samples to be used with the result from
# rejection_sample_parallel
tidy_samples_ind <- function(samples_list) {

  nc <- length(samples_list)
  rr <- nrow(samples_list[[1]])
  cc <- ncol(samples_list[[1]])

  # turn into array of samples
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- samples_list[[c]]
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                                 "intercept" = "...1",
                                 "vfi" = "...2",
                                 "tfi" = "...3",
                                 "slope" = "...4",
                                 "wind" = "...5",
                                 "steps" = "...6")

  return(draws_arr2)
}

# Constants --------------------------------------------------------------

# dir to load files
target_dir <- file.path("files", "pseudolikelihoods")

results_names <- list.files(target_dir)
samples_names <- results_names[grep("posterior_samples", results_names)]
strsplit(samples_names[1], split = "-posterior_samples.rds")[[1]]

rrr <- readRDS(file.path(target_dir, samples_names[1]))
str(rrr)
rrr <- do.call("rbind", rrr)
rrr <- tidy_samples_ind(rrr)
str(rrr)

%>% tidy_samples_ind()
str(rrr)




