# This script can be used to test how well our simulator works in a simple scenario
# where we have two uncorrelated confounders, with linear, time-invariant effects.
# We run exactly 10 models and save each results data frame separately.
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
k <- 2

# number of observed waves
T_waves <- 5

# number of initial waves to discard before keeping the observed panel
burn_in <- 20

# covariance of the observed base confounders
Omega11 <- diag(2)

# sample baseline Delta matrix at t = 1
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
  1.00, 0.30,
  0.30, 1.00
), nrow = 2, byrow = TRUE)

# create a dataframe for checking
out <- simulate_panel_data(
  N = 5000,
  T = T_waves,
  Phi = Phi,
  Omega11 = Omega11,
  Sigma = Sigma,
  Delta_list = Delta_list,
  burn_in = burn_in,
  eig_tol = 1e-10
)

out

# ------------------------------- shared analysis choices -------------------------------

exclude <- NULL
free_loadings <- FALSE
bootstrap_B <- 25
cores <- 4

# ------------------------------- XGB arguments -------------------------------

tune_xgb <- TRUE
xgb_tuning <- NULL

xgb_tune_args <- list(
  tuning_grid = expand.grid(
    eta = c(0.01, 0.03, 0.05, 0.10),
    max_depth = c(2, 3, 4, 6),
    min_child_weight = c(1, 3, 5, 7),
    subsample = c(0.7, 0.8, 1.0),
    colsample_bytree = c(0.7, 0.8, 1.0)
  ),
  cv_folds = 5,
  nrounds_max = 600,
  early_stopping_rounds = 30,
  nthread = 1,
  seed = 123
)

residualizer_args_xgb <- list(
  oof_folds = 5,
  nthread = 1,
  seed = 123
)

# linear residualizer should use the same OOF setup as XGB
residualizer_args_linear <- list(
  oof_folds = 5,
  seed = 123
)

# for residualizer = "none"
residualizer_args_none <- list()

# -------------------- model 1: CLPM without confounders --------------------

results_none_clpm_noconf <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "none",
  sem_model = "clpm",
  confounder_order = 0,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_none
)

df_none_clpm_noconf <- results_none_clpm_noconf$results

saveRDS(
  df_none_clpm_noconf,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "01_none_CLPM_no_confounders_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 2: CLPM with confounders in outcome model --------------------

results_none_clpm_withconf <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "none",
  sem_model = "clpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_none
)

df_none_clpm_withconf <- results_none_clpm_withconf$results

saveRDS(
  df_none_clpm_withconf,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_none_CLPM_with_confounders_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 3: linear residualizer + CLPM --------------------

results_linear_clpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "linear",
  sem_model = "clpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_linear
)

df_linear_clpm <- results_linear_clpm$results

saveRDS(
  df_linear_clpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "03_linearresid_CLPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 4: linear residualizer + RI-CLPM --------------------

results_linear_riclpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "linear",
  sem_model = "riclpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_linear
)

df_linear_riclpm <- results_linear_riclpm$results

saveRDS(
  df_linear_riclpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "04_linearresid_RICLPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 5: linear residualizer + DPM --------------------

results_linear_dpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "linear",
  sem_model = "dpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_linear
)

df_linear_dpm <- results_linear_dpm$results

saveRDS(
  df_linear_dpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "05_linearresid_DPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 6: xgb residualizer + CLPM --------------------

results_xgb_clpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "xgb",
  sem_model = "clpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = tune_xgb,
  xgb_tuning = xgb_tuning,
  xgb_tune_args = xgb_tune_args,
  residualizer_args = residualizer_args_xgb
)

df_xgb_clpm <- results_xgb_clpm$results

saveRDS(
  df_xgb_clpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "06_xgbresid_CLPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 7: xgb residualizer + RI-CLPM --------------------

results_xgb_riclpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "xgb",
  sem_model = "riclpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = tune_xgb,
  xgb_tuning = xgb_tuning,
  xgb_tune_args = xgb_tune_args,
  residualizer_args = residualizer_args_xgb
)

df_xgb_riclpm <- results_xgb_riclpm$results

saveRDS(
  df_xgb_riclpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "07_xgbresid_RICLPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 8: xgb residualizer + DPM --------------------

results_xgb_dpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "xgb",
  sem_model = "dpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = tune_xgb,
  xgb_tuning = xgb_tuning,
  xgb_tune_args = xgb_tune_args,
  residualizer_args = residualizer_args_xgb
)

df_xgb_dpm <- results_xgb_dpm$results

saveRDS(
  df_xgb_dpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "08_xgbresid_DPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 9: RI-CLPM only --------------------

results_none_riclpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "none",
  sem_model = "riclpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_none
)

df_none_riclpm <- results_none_riclpm$results

saveRDS(
  df_none_riclpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "09_none_RICLPM_constant_2c_linear_5000.rds"
  )
)

# -------------------- model 10: DPM only --------------------

results_none_dpm <- run_simulation_study(
  reps = 100,
  N = 5000,
  T = T_waves,
  burn_in = burn_in,
  k = k,
  Phi = Phi,
  Sigma = Sigma,
  Omega11 = Omega11,
  Delta_list = Delta_list,
  residualizer = "none",
  sem_model = "dpm",
  confounder_order = 1,
  exclude = exclude,
  free_loadings = free_loadings,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = 2024,
  cores = cores,
  base_seed = 1234,
  tune_xgb = FALSE,
  xgb_tuning = NULL,
  xgb_tune_args = list(),
  residualizer_args = residualizer_args_none
)

df_none_dpm <- results_none_dpm$results

saveRDS(
  df_none_dpm,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "10_none_DPM_constant_2c_linear_5000.rds"
  )
)

# optional: collect all data frames in one object as well
all_model_dfs <- list(
  df_none_clpm_noconf = df_none_clpm_noconf,
  df_none_clpm_withconf = df_none_clpm_withconf,
  df_linear_clpm = df_linear_clpm,
  df_linear_riclpm = df_linear_riclpm,
  df_linear_dpm = df_linear_dpm,
  df_xgb_clpm = df_xgb_clpm,
  df_xgb_riclpm = df_xgb_riclpm,
  df_xgb_dpm = df_xgb_dpm,
  df_none_riclpm = df_none_riclpm,
  df_none_dpm = df_none_dpm
)

saveRDS(
  all_model_dfs,
  file = here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "00_all_10_model_dataframes_constant_2c_linear_5000.rds"
  )
)