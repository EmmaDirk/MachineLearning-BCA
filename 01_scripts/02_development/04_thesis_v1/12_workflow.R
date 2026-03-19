# This master script shows the simplified final interface:
# - one data-generating Delta trajectory
# - one chosen residualiser
# - one chosen SEM
# - one output list with the long results frame plus the study objects
#
# This version is the parallel master script:
# - it is intended to be the script you run for larger studies
# - it uses the parallel branch in run_simulation_study()
# - it keeps XGB-related objects empty when the chosen residualiser is not XGB
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
set.seed(1234)

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
residualizer <- "linear"   # "none", "linear", or "xgb"
sem_model    <- "clpm"     # "clpm", "riclpm", or "dpm"

# analyst-side confounder choices shared across the pipeline where relevant
confounder_order <- 1
exclude <- NULL

# only relevant for RI-CLPM and DPM
free_loadings <- FALSE

# number of bootstrap resamples used for two-stage methods
bootstrap_B <- 50

# choose how many worker processes to use
cores <- 6

# ------------------------------- XGB-related arguments -------------------------------

# when residualizer is not XGB, leave XGB-related objects empty
# the simulation function now checks this explicitly
tune_xgb <- FALSE
xgb_tuning <- NULL
xgb_tune_args <- list()

# extra arguments passed only to the chosen residualiser during fitting
# for the linear residualiser, no extra arguments are needed here
residualizer_args <- list()

# ------------------------------- run the study -------------------------------

results_sim <- run_simulation_study(
  reps = 200,
  N = 5000,
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
save(results_sim, file = here("03_output", "02_thesis", "01_tests", "01_data","BCALM_CLPM_constant_1c_linear.RData"))

str(results_sim)
