# This master script shows the simplified final interface:
# - one data-generating Delta trajectory
# - one chosen residualiser
# - one chosen SEM
# - one output list with the long results frame plus the study objects
#
# This version is the parallel master script:
# - it is intended to be the script you run for larger studies
# - it uses the parallel branch in run_simulation_study()
# - here we use XGB residualisation and fit a CLPM
# -------------------------------------------------------------------------------------------------

library(here)

# load the study scripts
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "01_delta_sampler.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "02_delta_trajectory.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "04_lavaan_model_string_builder.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "05_residualisers.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "06_model_fitters.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "07_bootstrap_helpers.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "08_fit_stat_extractors.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "09_one_replication_wrapper.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "10_simulation_function.R"))

# reproducibility
set.seed(1233)

# ------------------------------- data-generating mechanism -------------------------------

# number of observed base confounders
k <- 1

# number of waves
T_waves <- 5

# covariance of the observed base confounders
Omega11 <- diag(1)

# sample baseline Delta matrix at t = 1
Delta1 <- sample_delta_1(
  k = k,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

# choose one Delta trajectory for this study
Delta_list <- generate_Delta_constant(
  Delta1 = Delta1$Delta,
  T = T_waves
)

# lag matrix used in the data-generating mechanism
Phi <- matrix(c(
  0.20, 0.00,
  0.10, 0.20
), nrow = 2, byrow = TRUE)

# target observed covariance of (X_t, Y_t)
Sigma <- matrix(c(
  1.00, 0.10,
  0.10, 1.00
), nrow = 2, byrow = TRUE)

# ------------------------------- analysis choices -------------------------------

# choose exactly one residualiser and exactly one SEM
residualizer <- "xgb"      # "none", "linear", or "xgb"
sem_model    <- "clpm"     # "clpm", "riclpm", or "dpm"

# analyst-side confounder choices shared across the pipeline where relevant
confounder_order <- 1
exclude <- NULL

# relevant for RI-CLPM and DPM
free_loadings <- FALSE

# bootstrap SEs are only relevant for two-stage methods
bootstrap_B <- 25

# choose how many worker processes to use
cores <- 7

# ------------------------------- XGB arguments -------------------------------

# cheap tuning setup
tune_xgb <- TRUE
xgb_tuning <- NULL

xgb_tune_args <- list(
  tuning_grid = expand.grid(
    eta = c(0.05, 0.10),
    max_depth = c(2, 3),
    min_child_weight = c(1, 5),
    subsample = 0.8,
    colsample_bytree = 0.8
  ),
  cv_folds = 3,
  nrounds_max = 150,
  early_stopping_rounds = 10,
  nthread = 1,
  seed = 123
)

# cheap fitting / residualising setup
residualizer_args <- list(
  oof_folds = 2,
  nthread = 1,
  seed = 123
)

# ------------------------------- run the study -------------------------------

results_sim <- run_simulation_study(
  reps = 100,
  N = 300,
  T = T_waves,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = residualizer,
  sem_model = sem_model,
  confounder_order = confounder_order,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = tune_xgb,
  xgb_tuning = xgb_tuning,
  xgb_tune_args = xgb_tune_args,
  residualizer_args = residualizer_args
)

# save the results
saveRDS(
  results_sim,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "03_BCAXGB_CLPM_constant_1c_linear_300.rds"
  )
)