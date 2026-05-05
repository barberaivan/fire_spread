# Fire Spread Modelling — Codebase Inventory

> Purpose: Feed this document to an LLM to design the architecture of a new,
> clean repo. Written from static analysis (file headers, grep patterns, no
> execution). Date: 2026-05-05.

---

## 1. Project Overview

**Domain:** Fire regime modelling for north-western Patagonia (Nahuel Huapi
National Park area). Part of a PhD thesis being converted into a journal paper,
with a future production side.

**Three fire sub-models:**

| Model | Purpose | Status |
|-------|---------|--------|
| Ignition | Probability of fire ignition per unit area per unit time | Fitted, canonical version in `ignition-escape_FWIZ/` |
| Escape | Probability a fire escapes to become a large fire | Fitted alongside ignition |
| Spread | Fire spread simulation, pixel-level burn probability | Most complex; fitted in `spread/` |

**Language stack:**
- R — all model fitting, simulation, data processing
- C++ via `Rcpp` — spread simulation engine (`FireSpread` package, external repo)
- One `.cpp` file inside this repo: `spread/sample_triplets_weighted.cpp`

**Working directory convention:** All scripts assume `fire_spread/` as the
working directory. Paths like `file.path("spread", "...")` are relative to that
root.

**Key external dependency:** `../FireSpread` — a separate R package wrapping
C++ spread functions. Referenced as:
```r
source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
library(FireSpread)
```
This path (`tests/testthat/`) is fragile — an important tech-debt item.

**Science vs. production distinction:**
- Fitted spread model parameters → constants usable by the production simulator
- Fire regime simulations → scientific output (many runs, analysis of results)
- Ignition/escape model outputs → also constants for production

---

## 2. Directory Map

```
fire_spread/
├── spread/                     # Spread model fitting — largest module (~26 R + 1 cpp)
├── ignition-escape_FWIZ/       # Ignition & escape model fitting — CANONICAL version
├── ignition-escape/            # OLD — superseded by _FWIZ
├── ignition-escape_FWIZ_OLD-EXPONENTIAL/  # OLD — superseded
├── fire regime simulations/    # Integration of all 3 models; scientific simulation runs
├── flammability indices/       # Preprocessing: NDVI/elevation → VFI/TFI indices
├── weather/                    # FWI and wind preprocessing; fortnight aggregation
│   └── FWI projections/        # CMIP6 projected FWI processing (2050, 2090)
├── data/                       # Raw inputs (gitignored heavy files)
│   ├── focal fires data/       # Landscape rasters + GEE exports for each focal fire
│   │   └── landscapes/         # Prepared landscape arrays (output of landscapes_preparation.R)
│   ├── fwi_daily_1998-2022/    # Daily FWI rasters (24 km and 52 km resolution)
│   ├── fwi_projections/        # CMIP6 future FWI data
│   ├── flammability indices/   # Fitted flammability index parameters (.rds)
│   ├── ignition/               # Ignition point data
│   ├── pnnh_images/            # PNNH landscape rasters (6-layer; wind, FWI, VFI, TFI)
│   ├── protected_areas/        # Shapefiles: park boundaries, buffers
│   ├── simulations/            # Raw simulation outputs
│   └── wind_2020-2022/         # Hourly wind data
├── files/                      # Computed outputs (gitignored)
│   ├── posterior_samples_stage1_smc/   # Stage-1 fire-wise SMC posteriors
│   ├── hierarchical_model_FWIZ_SMC/    # Fitted hierarchical spread model
│   ├── ignition_FWIZ/                  # Fitted ignition/escape model samples
│   ├── fire_regime_simulation_FWIZ/    # Simulation runs output
│   ├── landscape_flammability/         # Precomputed flammability maps
│   └── ... (many intermediate output dirs, mostly derivable by re-running)
└── dump/                       # Archived / exploratory code — do not carry forward
```

---

## 3. Canonical Scripts

The "one true version" of each role. Everything else is legacy or exploratory.

### 3a. Function Libraries (Layer 1 — no upstream script dependencies)

| File | Role | What it exports |
|------|------|-----------------|
| `flammability indices/flammability_indices_functions.R` | VFI/TFI computation + NDVI detrending | `ndvi_detrend()`, `vfi()`, `tfi()` functions; loads `data/flammability indices/flammability_indices.rds` and `ndvi_detrender_model.rds` |
| `weather/fortnight_functions.R` | 14-day fortnight indexing | `date2fort()` and reference table; origin fixed at 1996 for FWI compatibility |
| `spread/mcmc_functions.R` | Core MCMC utilities for hierarchical spread model | Gibbs updates, inverse-Wishart, MH for random effects, scaled-logit-normal steps |
| `spread/mcmc_functions_SMC.R` | Extends `mcmc_functions.R` for SMC stage-1 | Modified `update_ranef()` using logistic(0,1) prior from SMC samples |
| `../FireSpread/tests/testthat/R_spread_functions.R` | R wrappers for C++ spread engine | `rast_from_mat()`, constants; **external repo** |

### 3b. Data Preprocessing (Layer 2a — produce intermediate data)

| File | Role | Reads | Writes |
|------|------|-------|--------|
| `flammability indices/flammability_indices.R` | Fits VFI/TFI models via Stan logistic regression | Raw NDVI rasters, fire data from `data/` | `data/flammability indices/flammability_indices.rds`, `ndvi_detrender_model.rds` |
| `weather/FWI standardize and aggregate by fortnight.R` | Detrends FWI to temporal anomalies; aggregates by fortnight | `data/fwi_daily_1998-2022/24km/*.tif` | `data/fwi_fortnights_*_standardized.tif` |
| `weather/FWI fortnight matrix for spread and lengthscale estimation.R` | Inherits from `FWI temporal scale.R`; uses FWI **anomalies** instead of raw values; interpolates FWI at ignition points | FWI standardized rasters + ignition shapefiles | Lagged FWI matrix CSV for model fitting |
| `weather/FWI projections/standardize and aggregate projections by fortnight (2050 and 2090).R` | Processes CMIP6 projected FWI using modern-period calibration | `data/fwi_projections/` | `data/fwi_projections/fwi_fortnights_standardized/` |

### 3c. Model Fitting (Layer 2b — produce fitted model objects)

| File | Role | Reads | Writes |
|------|------|-------|--------|
| `spread/landscapes_preparation.R` | Creates 6-layer landscape arrays (VFI, TFI, elev, wind dir, wind speed, FWI) for each focal fire | `data/focal fires data/raw data from GEE/`, `data/pnnh_images/`, windninja outputs; sources `R_spread_functions.R` + `flammability_indices_functions.R` | `data/focal fires data/landscapes/*.rds` |
| `spread/sampling_fire_wise_posteriors_(stage1)_SMC.R` | Stage-1: samples fire-wise posterior distributions via ABC-SMC (Del Moral et al. 2011); uses hard ABC kernel | Landscape arrays from `data/focal fires data/landscapes/`, `flammability_indices.rds`; sources `R_spread_functions.R`, compiles `sample_triplets_weighted.cpp` | `files/posterior_samples_stage1_smc/*.rds` |
| `spread/hierarchical model fitting_FWIZ2_SMC.R` | Stage-2: fits hierarchical Bayesian spread model via custom MCMC (Gibbs + MH); uses SMC stage-1 samples as proposals | `files/posterior_samples_stage1_smc/`, FWI fortnight matrix; sources `mcmc_functions_SMC.R` + `flammability_indices_functions.R` | `files/hierarchical_model_FWIZ_SMC/*.rds` |
| `ignition-escape_FWIZ/ignition-escape_analyses.R` | Fits ignition (negative binomial) and escape (logistic) models via Stan | `data/ignition/ignition_size_data.csv`, ignition shapefiles, FWI data; sources `flammability_indices_functions.R` + `fortnight_functions.R` | `files/ignition_FWIZ/ignition_model_samples.rds`, `escape_model_samples.rds` |

### 3d. C++ component

| File | Role |
|------|------|
| `spread/sample_triplets_weighted.cpp` | Compiled via `sourceCpp()`; samples triplets of (burned area, fire shape similarity, steps) weighted by likelihood — core of ABC-SMC stage-1 |

### 3e. Integration & Simulation (Layer 3)

| File | Role | Reads | Writes |
|------|------|-------|--------|
| `fire regime simulations/fire_regime_simulations_FWIZ.R` | Full fire regime simulator: integrates ignition, escape, spread models; runs many simulations for scientific analysis; recalibrates spread parameters for PNNH | `files/hierarchical_model_FWIZ_SMC/`, `files/ignition_FWIZ/`, `data/pnnh_images/pnnh_*.tif`, FWI fortnight rasters; sources `R_spread_functions.R` + `flammability_indices_functions.R` + `fortnight_functions.R` | `files/fire_regime_simulation_FWIZ/*.rds` |
| `fire regime simulations/fire_probability_maps_single-models_FWIZ.R` | Generates static fire probability maps from single model runs | `files/fire_regime_simulation_FWIZ/` | Figures |
| `fire regime simulations/plots.R` | Visualization utilities for regime simulation results | `files/fire_regime_simulation_FWIZ/` | Figures (A4 format) |

---

## 4. Full Canonical Pipeline

```
Raw data (GEE exports, FWI tifs, fire shapefiles)
    │
    ├─ flammability indices/flammability_indices.R
    │       └─ data/flammability indices/flammability_indices.rds
    │              ndvi_detrender_model.rds
    │
    ├─ weather/FWI standardize and aggregate by fortnight.R
    │       └─ data/fwi_daily_1998-2022/24km/fwi_fortnights_*_standardized.tif
    │
    ├─ spread/landscapes_preparation.R
    │       └─ data/focal fires data/landscapes/*.rds  (one per fire)
    │
    ├─ ignition-escape_FWIZ/ignition-escape_analyses.R
    │       └─ files/ignition_FWIZ/{ignition,escape}_model_samples.rds
    │
    ├─ spread/sampling_fire_wise_posteriors_(stage1)_SMC.R
    │       └─ files/posterior_samples_stage1_smc/*.rds
    │
    ├─ spread/hierarchical model fitting_FWIZ2_SMC.R
    │       └─ files/hierarchical_model_FWIZ_SMC/*.rds
    │
    └─ fire regime simulations/fire_regime_simulations_FWIZ.R
            └─ files/fire_regime_simulation_FWIZ/*.rds
                    └─ plotting scripts → figures/
```

---

## 5. Dependency Graph (source/import connections)

### Function library sourcing

```
flammability_indices_functions.R  ←── flammability_indices.R (fits and saves params)
    ↑ sourced by:
    spread/landscapes_preparation.R
    spread/hierarchical model fitting_FWIZ2_SMC.R
    ignition-escape_FWIZ/ignition-escape_analyses.R
    fire regime simulations/fire_regime_simulations_FWIZ.R
    fire regime simulations/fire_probability_maps_single-models_FWIZ.R
    fire regime simulations/plots_defensa_balconGut.R

weather/fortnight_functions.R
    ↑ sourced by:
    weather/FWI standardize and aggregate by fortnight.R
    weather/FWI fortnight matrix for spread and lengthscale estimation.R
    weather/FWI projections/standardize and aggregate projections by fortnight (2050 and 2090).R
    ignition-escape_FWIZ/ignition-escape_analyses.R
    ignition-escape_FWIZ/ignition-escape_analyses-with-sizeclass.R
    ignition-escape_FWIZ/ignition-escape_analyses-trying temporal models.R
    fire regime simulations/fire_regime_simulations_FWIZ.R

spread/mcmc_functions.R
    ↑ sourced by:
    spread/hierarchical model fitting.R  (legacy)
    spread/hierarchical model fitting_FWIZ.R  (legacy)
    spread/hierarchical model fitting_FWIZ2.R  (legacy)
    spread/hierarchical model fitting_FWIZ2_extra steps models.R

spread/mcmc_functions_SMC.R  (sources mcmc_functions.R)
    ↑ sourced by:
    spread/hierarchical model fitting_FWIZ2_SMC.R  ← CANONICAL

../FireSpread/tests/testthat/R_spread_functions.R
    ↑ sourced by:
    spread/landscapes_preparation.R
    spread/sampling_fire_wise_posteriors_(stage1)_SMC.R
    spread/sampling_fire_wise_posteriors_(stage1).R  (legacy)
    spread/hierarchical model fitting_FWIZ2.R  (legacy)
    spread/hierarchical model fitting_FWIZ2_SMC.R
    spread/sampling_fire_wise_posteriors_SURROGATES.R  (legacy)
    fire regime simulations/fire_regime_simulations_FWIZ.R
    fire regime simulations/plots_defensa_balconGut.R

spread/sample_triplets_weighted.cpp
    ↑ compiled (sourceCpp) by:
    spread/sampling_fire_wise_posteriors_(stage1)_SMC.R  ← CANONICAL
```

### Key data flow (readRDS / rast / vect / read.csv)

```
data/flammability indices/flammability_indices.rds
    → loaded inside flammability_indices_functions.R (at source time)

data/flammability indices/ndvi_detrender_model.rds
    → loaded inside flammability_indices_functions.R (at source time)

data/fwi_daily_1998-2022/24km/fwi_daily_19970701_20230630.tif  (rast)
    → weather/FWI standardize and aggregate by fortnight.R

data/focal fires data/landscapes/*.rds
    → spread/sampling_fire_wise_posteriors_(stage1)_SMC.R

files/posterior_samples_stage1_smc/*.rds
    → spread/hierarchical model fitting_FWIZ2_SMC.R

files/hierarchical_model_FWIZ_SMC/*.rds
    → fire regime simulations/fire_regime_simulations_FWIZ.R

files/ignition_FWIZ/ignition_model_samples.rds
files/ignition_FWIZ/escape_model_samples.rds
    → fire regime simulations/fire_regime_simulations_FWIZ.R

data/pnnh_images/pnnh_data_30m_buff_10000_std.tif  (rast, 6-layer)
    → fire regime simulations/fire_regime_simulations_FWIZ.R

data/protected_areas/apn_limites.shp  (vect — most referenced shapefile)
    → spread/landscapes_preparation.R
    → ignition-escape_FWIZ/ignition-escape_analyses.R

data/ignition_points_checked_with_date-fort-fwiz2.shp  (vect — canonical ignition points)
    → ignition-escape_FWIZ/ignition-escape_analyses.R

data/climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv  (canonical climate-per-fire CSV)
    → spread/hierarchical model fitting_FWIZ2_SMC.R
    → spread/hierarchical model fitting_FWIZ2.R

data/ignition/ignition_size_data.csv
    → ignition-escape_FWIZ/ignition-escape_analyses.R

External shapefile (outside repo):
/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp
    → data preparation / ignition scripts
```

---

## 6. Legacy and Deprecated Files

These should **not** be carried into the new repo (or archived in a `legacy/` folder):

### Superseded by FWIZ variants
- `ignition-escape/ignition-escape_analyses.R` — old (pre-FWIZ) version
- `ignition-escape_FWIZ_OLD-EXPONENTIAL/ignition-escape_analyses.R` — old parameterization
- `fire regime simulations/fire_regime_simulations.R` — non-FWIZ version
- `fire regime simulations/fire_probability_maps_single-models.R` — non-FWIZ

### Superseded by FWIZ2 / SMC variants in spread/
- `spread/hierarchical model fitting.R` — base version
- `spread/hierarchical model fitting_FWIZ.R` — intermediate
- `spread/hierarchical model fitting_FWIZ2.R` — dissertation; SMC version is canonical
- `spread/sampling_fire_wise_posteriors_(stage1).R` — dissertation (rejection sampling)

### Development / exploration artifacts
- `spread/sampling_fire_wise_posteriors_IMPORTANCE.R` — importance sampling approach, not used
- `spread/sampling_fire_wise_posteriors_SURROGATES.R` — surrogate likelihood approach, not used
- `spread/smc_0.R` — development of SMC algorithm, merged into canonical scripts
- `spread/mcmc_functions_BACKUP-OLD.R` — backup
- `spread/mcmc_functions_test.R` — test/dev
- `spread/hierarchical model fitting_FWIZ2_extra steps models.R` — model exploration
- `spread/exploring fire-wise posteriors to define model structure.R` — exploration
- `spread/exploring fire-wise posteriors to define model structure_FWIZ.R` — exploration
- `spread/gibbs linear regression.R` — reference implementation
- `spread/gibbs_testing.R` — development notes
- `spread/FireSpread package interactive testing.R` — interactive testing
- `spread/jacobian.R`, `spread/log-uniform density.R`, `spread/scaled-logit-normal density.R`, `spread/normal_pseudolikelihoods.R` — density / math scratch scripts
- `spread/sampling_fire_wise_posteriors_IMPORTANCE.R` — not used

### Old weather scripts
- `weather/fortnight_functions_OLD.R` — superseded by `fortnight_functions.R`
- `weather/climate temporal scale.R` — explored temp/precip before choosing FWI
- `weather/FWI temporal scale.R` — exploration; superseded by the fortnight matrix script
- `weather/wind temporal scale.R` — exploration
- `weather/FWI projections/standardize and aggregate projections by fortnight (2050 and 2090)_OLD.R`
- `weather/FWI projections/explore available data_old.R`

### Old flammability scripts
- `flammability indices/comparing ndvi_max and ndvi_dyn.R` — exploration
- `flammability indices/northing importance function.R` — exploration

### Data helper scripts (post-processing in wrong location)
- `data/simulations/burned_mat.R`, `burned_mat_all.R`, `times.R`, `times_all.R`
  — post-process raw simulation output; belong in `fire regime simulations/`

### Dump folder
All files in `dump/` — archived exploration, not used:
- 3 versions of hierarchical model fitting (`_01`, `_02`, `_03`)
- 2 versions of landscapes_preparation
- 2 versions of FireSpread interactive testing
- exploratory posterior and overlap scripts
- `notes on estimation and roadmap.txt` — useful reading for context

---

## 7. R Packages Used

Collected from `library()` calls across canonical scripts. A new repo should
document these explicitly.

**Core data/spatial:**
`terra`, `tidyverse`, `lubridate`, `dplyr`, `stringr`

**Bayesian / statistical:**
`rstan`, `posterior`, `bayesplot`, `bayestestR`, `tidybayes`, `LaplacesDemon`,
`invgamma`, `MASS`, `truncnorm`, `truncreg`, `logitnorm`, `mgcv`, `DHARMa`,
`HDInterval`, `loo`, `brms`, `abind`

**MCMC / sampling:**
`randtoolbox` (Sobol sequences), `trialr` (LKJ correlation)

**Parallelization:**
`foreach`, `doMC`

**Spread engine:**
`FireSpread` (external), `Rcpp`

**Visualization:**
`viridis`, `deeptime`, `patchwork`, `tidyterra`, `ggplot2` (via tidyverse)

**Utilities:**
`microbenchmark`, `rvinecopulib` (in legacy surrogate scripts)

---

## 8. Key Data Files to Track

Files that are gitignored (heavy) but required to reproduce results:

### Input rasters (`data/`)
| File | Description |
|------|-------------|
| `data/fwi_daily_1998-2022/24km/fwi_daily_19970701_20230630.tif` | Daily FWI, 25 years, 24 km resolution |
| `data/fwi_daily_1998-2022/24km/fwi_fortnights_*_standardized.tif` | Processed fortnight FWI anomalies |
| `data/pnnh_images/pnnh_data_30m_buff_10000_std.tif` | PNNH 6-layer landscape raster (30 m) |
| `data/pnnh_images/wind_direction_*.tif`, `wind_speed_*.tif` | Wind layers |
| `data/fwi_projections/fwi_fortnights_standardized/` | CMIP6 2050/2090 FWI projections |
| `data/focal fires data/landscapes/*.rds` | Per-fire landscape arrays (output of `landscapes_preparation.R`) |

### Input vectors (`data/`)
| File | Description |
|------|-------------|
| `data/protected_areas/apn_limites.shp` | National park boundary |
| `data/protected_areas/pnnh_buff_10000.shp` | PNNH + 10 km buffer |
| `data/ignition_points_checked_with_date-fort-fwiz2.shp` | Canonical ignition points with fortnight dates |
| `data/patagonian_fires_spread.shp` | Fire perimeters used for spread fitting |
| `/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp` | **External** — all Patagonian fires |

### Input CSVs (`data/`)
| File | Description |
|------|-------------|
| `data/climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv` | Canonical FWI per fire |
| `data/ignition/ignition_size_data.csv` | Ignition/escape training data |

### Key computed outputs (`files/` — derivable by running pipeline)
| File | Description |
|------|-------------|
| `files/hierarchical_model_FWIZ_SMC/*.rds` | Fitted spread model — **production constant** |
| `files/ignition_FWIZ/ignition_model_samples.rds` | Fitted ignition model — **production constant** |
| `files/ignition_FWIZ/escape_model_samples.rds` | Fitted escape model — **production constant** |
| `files/posterior_samples_stage1_smc/*.rds` | Fire-wise posteriors (stage-1 intermediate) |

---

## 9. Tech Debt / Refactoring Notes for New Repo

1. **`spread/landscapes_preparation.R` is a loop, not a function.** The CLAUDE.md
   explicitly flags this: it should become a function that can create any
   landscape, not a hard-coded loop over focal fires.

2. **`spread/hierarchical model fitting_FWIZ2_SMC.R` is monolithic (~2800 lines).**
   Contains both MCMC algorithm logic and data manipulation. The MCMC core is
   already extracted into `mcmc_functions_SMC.R`, but inline sections in the
   fitting script should move to functions.

3. **External `R_spread_functions.R` sourced from `tests/testthat/` path.**
   This file lives inside the `FireSpread` package test suite — a fragile
   location. In a new repo, this should be a proper exported function or a
   vendored utility file.

4. **`data/simulations/*.R` scripts live in `data/` but belong in
   `fire regime simulations/`.** They post-process raw simulation arrays and
   are conceptually part of the simulation pipeline.

5. **Multiple parallel naming conventions.** Scripts with `_FWIZ`, `_FWIZ2`,
   `_SMC` suffixes exist because the project evolved. The new repo should use
   a single version with version history in git.

6. **Working directory assumption is implicit.** All `file.path(...)` calls
   assume `fire_spread/` as CWD. A new repo should make this explicit via a
   config file or `.Rprofile`.

7. **`windninja_cli_fire_spread_files` path is hardcoded** at
   `/home/ivan/windninja_cli_fire_spread_files` in `landscapes_preparation.R`.
   This is a machine-specific absolute path.

---

## 10. What the New Repo Should Contain

Based on canonical scripts only:

```
new_repo/
├── R/
│   ├── flammability_indices_functions.R    # function library
│   ├── fortnight_functions.R               # function library
│   ├── mcmc_functions.R                    # MCMC utilities
│   └── mcmc_functions_smc.R               # SMC variant
├── src/
│   └── sample_triplets_weighted.cpp        # only C++ file
├── data_prep/
│   ├── flammability_indices.R              # fit VFI/TFI
│   ├── fwi_standardize.R                  # FWI fortnight standardization
│   ├── fwi_fortnight_matrix.R             # lagged FWI matrix
│   ├── fwi_projections.R                  # CMIP6 processing
│   └── landscapes_preparation.R           # landscape arrays → function
├── spread/
│   ├── stage1_smc.R                       # fire-wise posteriors
│   └── hierarchical_fit.R                 # hierarchical MCMC
├── ignition_escape/
│   └── fit.R                              # ignition + escape models
├── fire_regime/
│   ├── simulate.R                         # main simulation
│   ├── probability_maps.R
│   └── plots.R
├── data/                                  # (gitignored heavy files)
├── files/                                 # (gitignored outputs)
└── CLAUDE.md
```
