# =========================================================================================
# smoke_test_model_set.R
# Purpose:
#   Very cheap mechanical test of the full shared-data simulation workflow.
#
# This script STILL tests:
# - burn-in
# - bootstrap SEs
# - XGBoost tuning
# - XGBoost residualisation
# - the shared model-set runner
#
# But it is made deliberately cheap by using:
# - very few replications
# - small N
# - small bootstrap B
# - tiny XGB tuning grid
# - few CV folds
# - few OOF folds
# - one saved combined dataframe only
# =========================================================================================

library(here)

# -----------------------------------------------------------------------------------------
# load the study scripts
# -----------------------------------------------------------------------------------------
source(here("01_scripts", "02_development", "05_thesis_v2", "00_packages.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "01_delta_sampler.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "02_delta_trajectory.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "04_lavaan_model_string_builder.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "05_residualisers.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "06_model_fitters.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "07_bootstrap_helpers.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "08_fit_stat_extractors.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "09_model_set_helpers.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "10_one_replication_wrapper.R"))
source(here("01_scripts", "02_development", "05_thesis_v2", "11_simulation_function.R"))

# -----------------------------------------------------------------------------------------
# reproducibility
# -----------------------------------------------------------------------------------------
set.seed(1233)

# -----------------------------------------------------------------------------------------
# data-generating mechanism
# -----------------------------------------------------------------------------------------
k <- 2
T_waves <- 5
burn_in <- 5                     # still tests burn-in, but cheap
Omega11 <- diag(2)

Delta1 <- sample_delta_1(
  k = k,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1,
  force_positive = TRUE
)

Delta_list <- generate_Delta_constant(
  Delta1 = Delta1$Delta,
  T = T_waves
)

Phi <- matrix(c(
  0.20, 0.00,
  0.10, 0.20
), nrow = 2, byrow = TRUE)

Sigma <- matrix(c(
  1.00, 0.30,
  0.30, 1.00
), nrow = 2, byrow = TRUE)

# optional quick sanity check of one generated dataset
sim_check <- simulate_panel_data(
  N = 40,
  T = T_waves,
  Phi = Phi,
  Omega11 = Omega11,
  Sigma = Sigma,
  Delta_list = Delta_list,
  burn_in = burn_in,
  eig_tol = 1e-10
)

print(names(sim_check))
print(head(sim_check$data))

# -----------------------------------------------------------------------------------------
# cheap smoke-test settings
# -----------------------------------------------------------------------------------------
reps <- 3                        # cheap but >1
N <- 60                          # small sample
bootstrap_B <- 4                 # still triggers bootstrap logic
cores <- 3                       # simplest and safest for a smoke test
exclude <- NULL
free_loadings <- FALSE

# -----------------------------------------------------------------------------------------
# cheap XGB settings
# still does real tuning + real XGB residualisation, just cheaply
# -----------------------------------------------------------------------------------------
tune_xgb <- TRUE
xgb_tuning <- NULL

xgb_tune_args <- list(
  tuning_grid = expand.grid(
    eta = c(0.10, 0.30),
    max_depth = c(2, 3),
    min_child_weight = 1,
    subsample = 0.8,
    colsample_bytree = 0.8
  ),
  cv_folds = 2,
  nrounds_max = 20,
  early_stopping_rounds = 5,
  nthread = 1,
  seed = 123
)

residualizer_args_xgb <- list(
  oof_folds = 2,
  nthread = 1,
  seed = 123
)

residualizer_args_linear <- list(
  oof_folds = 2,
  seed = 123
)

residualizer_args_none <- list()

# -----------------------------------------------------------------------------------------
# define a SMALL but COMPLETE model set
#
# This covers:
# - no residualisation
# - direct confounder-adjusted CLPM
# - linear residualisation
# - xgb residualisation
# - CLPM
# - RI-CLPM
# - DPM
#
# That is enough for a mechanical smoke test.
# -----------------------------------------------------------------------------------------
model_specs <- list(

  # 1) plain CLPM without confounders
  make_model_spec(
    name = "none_clpm_noconf",
    residualizer = "none",
    sem_model = "clpm",
    confounder_order = 0,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_none
  ),

  # 2) direct-adjusted CLPM
  make_model_spec(
    name = "none_clpm_withconf",
    residualizer = "none",
    sem_model = "clpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_none
  ),

  # 3) linear residualiser + CLPM
  make_model_spec(
    name = "linear_clpm",
    residualizer = "linear",
    sem_model = "clpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_linear
  ),

  # 4) linear residualiser + RI-CLPM
  make_model_spec(
    name = "linear_riclpm",
    residualizer = "linear",
    sem_model = "riclpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_linear
  ),

  # 5) xgb residualiser + CLPM
  make_model_spec(
    name = "xgb_clpm",
    residualizer = "xgb",
    sem_model = "clpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = tune_xgb,
    xgb_tuning = xgb_tuning,
    xgb_tune_args = xgb_tune_args,
    residualizer_args = residualizer_args_xgb
  ),

  # 6) xgb residualiser + DPM
  make_model_spec(
    name = "xgb_dpm",
    residualizer = "xgb",
    sem_model = "dpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = tune_xgb,
    xgb_tuning = xgb_tuning,
    xgb_tune_args = xgb_tune_args,
    residualizer_args = residualizer_args_xgb
  ),

  # 7) RI-CLPM only
  make_model_spec(
    name = "none_riclpm",
    residualizer = "none",
    sem_model = "riclpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_none
  ),

  # 8) DPM only
  make_model_spec(
    name = "none_dpm",
    residualizer = "none",
    sem_model = "dpm",
    confounder_order = 1,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    tune_xgb = FALSE,
    xgb_tuning = NULL,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args_none
  )
)

# -----------------------------------------------------------------------------------------
# run the smoke test
# -----------------------------------------------------------------------------------------
simulation_out <- run_simulation_model_set(
  reps = reps,
  N = N,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  model_specs = model_specs,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234
)

# -----------------------------------------------------------------------------------------
# combine all model outputs into ONE dataframe
# -----------------------------------------------------------------------------------------
all_model_df <- dplyr::bind_rows(
  simulation_out$results_by_model,
  .id = "model_name"
)

# inspect
print(dplyr::glimpse(all_model_df))
print(dplyr::count(all_model_df, model_name, model, residualizer))

# -----------------------------------------------------------------------------------------
# save ONE combined dataframe only
# -----------------------------------------------------------------------------------------
saveRDS(
  all_model_df,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "00_smoke_test_combined_model_dataframe.rds"
  )
)