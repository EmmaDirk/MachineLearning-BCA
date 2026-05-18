# =================================================================================================
#
# simulation runner for the final thesis scenarios
#
# Data-generating mechanism:
#
# - n_waves = 5 observed waves.
# - burn_in = 20 generated waves are discarded before analysis.
# - R2_total = 0.15 for the direct confounder contribution at the scaling wave.
# - R2_nl = 0.05 for the absolute interaction contribution at the scaling wave.
# - R2_new = 0.35 after the stepwise increase in direct confounder contribution.
# - step_at = 3, so the step begins at the third observed wave.
# - k = 5 baseline confounders.
# - Phi = {0.25, 0.00
#          0.10, 0.25}
# - Sigma = {1.00, 0.50
#            0.50, 1.00}
#
# Final thesis scenarios:
#
# 1. Constant linear confounding effects.
# 2. Constant interaction effects, with R2_nl = 0.05.
# 3. Stepwise rank-stable interaction effects, with R2_nl = 0.05,
#    R2_new = 0.35, and step_at = 3.
# 4. The same data-generating mechanism as scenario 3, but C1 and C2 are omitted
#    from every analyst-side adjustment.
#
# Sample sizes:
#
# - N = 2000.
# - N = 1000.
# - N = 300.
#
# Models fitted in every scenario:
#
# - CLPM without adjustment.
# - RI-CLPM without adjustment.
# - DPM without adjustment.
# - CLPM with observed linear confounder control.
# - RI-CLPM with observed linear confounder control.
# - DPM with observed linear confounder control.
# - CLPM after XGB residualisation.
# - RI-CLPM after XGB residualisation.
# - DPM after XGB residualisation.
#
# Output:
#
# - One .rds file is saved per scenario x sample-size condition.
# - One combined .rds file is saved after the full grid has finished.
#
# =================================================================================================

# ---- project paths -------------------------------------------------------------------------------

library(here)

# Set the root folder for this simulation study.
# Replace the text inside here::here() with the relative path to the study folder.

study_root <- here::here("[put your study root folder here]")

# Reusable directories.

script_dir <- file.path(study_root, "01_data_generation")
out_dir    <- file.path(study_root, "03_output")

# ---- logistics -----------------------------------------------------------------------------------

# Study script files.
# The master session and parallel workers source the same script versions.

simulation_worker_script_files <- c(
  "01_delta_sampler.R",
  "02_delta_trajectory.R",
  "03_simulate_panel_data.R",
  "04_lavaan_model_string_builder.R",
  "05_residualisers.R",
  "06_model_fitters.R",
  "07_bootstrap_helpers.R",
  "08_fit_stat_extractors.R",
  "09_model_set_helpers.R",
  "10_one_replication_wrapper.R",
  "11_simulation_function.R"
)

simulation_script_files <- c(
  "00_packages.R",
  simulation_worker_script_files
)

missing_script_files <- simulation_script_files[
  !file.exists(file.path(script_dir, simulation_script_files))
]

if (length(missing_script_files) > 0L) {
  stop(
    "Missing study scripts in ", script_dir, ": ",
    paste(missing_script_files, collapse = ", ")
  )
}

message("Study root: ", study_root)
message("Script directory: ", script_dir)
message("Output directory: ", out_dir)

# Source study scripts.
# All helper functions are loaded from the relative script directory.

for (f in simulation_script_files) {
  source(file.path(script_dir, f))
}

# Reproducibility.
# This seed controls deterministic setup outside the replication-level seeds.

set.seed(1233)

# Thread control.
# The study parallelises across replications, so each R process should use one
# BLAS/OpenMP thread to avoid thread oversubscription.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OMP_DYNAMIC = "FALSE",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  MKL_DYNAMIC = "FALSE",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
}

# ---- shared dgm ----------------------------------------------------------------------------------

# number of confounders
k <- 5

# number of observed waves
T_waves <- 5

# burnin
burn_in <- 20

# total variance explained by confounders at baseline
R2_total <- 0.15

# absolute nonlinear variance contribution for the nonlinear scenarios
R2_nl <- 0.05

# stepwise scenarios
R2_new <- 0.35
step_at <- 3L

# analyst-side omissions used only in scenario 4
omit_vars <- c("C1", "C2")

# varcov matrix of the main effects 
Omega11 <- matrix(c(
  1,    0,    0.05, 0.10, 0.15,
  0,    1,    0.20, 0.15, 0.10,
  0.05, 0.20, 1,    0.05, 0.10,
  0.10, 0.15, 0.05, 1,    0.15,
  0.15, 0.10, 0.10, 0.15, 1
), nrow = 5, byrow = TRUE)

# autoregressive and cross-lagged effects
Phi <- matrix(c(
  0.25, 0.00,
  0.10, 0.25
), nrow = 2, byrow = TRUE)

# varcov matrix of the observed X_t and Y_t
Sigma <- matrix(c(
  1.00, 0.50,
  0.50, 1.00
), nrow = 2, byrow = TRUE)

# ---- shared model specs --------------------------------------------------------------------------

# global 
reps <- 2000
bootstrap_B <- 40

# The runner is the only place where worker count is chosen.
# The default uses at most half of the available logical cores.
n_cores_detected <- parallel::detectCores(logical = TRUE)

if (is.na(n_cores_detected)) {
  cores <- 1L
} else {
  cores <- max(1L, floor(n_cores_detected / 2L))
}

# XGBoost and threaded math libraries are kept single-threaded inside each worker.
worker_nthread <- 1L

# XGB
tune_xgb <- TRUE
xgb_tuning <- NULL

xgb_tune_args <- list(
  tuning_grid = expand.grid(
    eta = c(0.02, 0.03, 0.05, 0.10),
    max_depth = c(1,2,3,4),
    min_child_weight = c(3, 5, 10),
    subsample = c(0.6, 0.8, 1.0),
    colsample_bytree = c(0.6, 0.8, 1.0)
  ),
  cv_folds = 8,
  nrounds_max = 500,
  early_stopping_rounds = 30,
  early_tune_wave = 2L,
  late_tune_wave = 4L,
  late_start_wave = 4L,
  nthread = worker_nthread,
  seed = 123
)

residualizer_args_xgb <- list(
  oof_folds = 8,
  use_oof = TRUE,
  nthread = worker_nthread,
  seed = 123
)

# elastic net 
# not used in the sim study, but arguments are required
tune_enet <- TRUE
enet_tuning <- NULL

enet_tune_args <- list(
  alpha_grid = seq(0, 1, by = 0.05),
  cv_folds = 5,
  seed = 123
)

residualizer_args_enet <- list(
  oof_folds = 5,
  seed = 123
)

# linear 
residualizer_args_linear <- list(
  oof_folds = 8,
  seed = 123
)

# none
residualizer_args_none <- list()

# ---- helper functions ----------------------------------------------------------------------------

# helper to build the model set for one scenario
# analyst_order describes what every adjustable model gets to observe
# exclude_general applies to every adjustable model
make_model_set <- function(
    analyst_order = 1L,
    exclude_general = NULL
) {
  
  list(
    # ---- no adjustment models --------------------------------------------------------------------
    make_model_spec(
      name = "clpm_no_adjustment",
      residualizer = "none",
      sem_model = "clpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    make_model_spec(
      name = "riclpm_no_adjustment",
      residualizer = "none",
      sem_model = "riclpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    make_model_spec(
      name = "dpm_no_adjustment",
      residualizer = "none",
      sem_model = "dpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    # ---- observed linear confounder control ------------------------------------------------------
    make_model_spec(
      name = "clpm_linear_confounders",
      residualizer = "none",
      sem_model = "clpm",
      sem_c_order = analyst_order,
      sem_exclude = exclude_general,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    make_model_spec(
      name = "riclpm_linear_confounders",
      residualizer = "none",
      sem_model = "riclpm",
      sem_c_order = analyst_order,
      sem_exclude = exclude_general,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    make_model_spec(
      name = "dpm_linear_confounders",
      residualizer = "none",
      sem_model = "dpm",
      sem_c_order = analyst_order,
      sem_exclude = exclude_general,
      residualizer_c_order = 0,
      residualizer_exclude = NULL,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = NULL,
      enet_tuning = NULL,
      tune_xgb = FALSE,
      tune_enet = FALSE,
      xgb_tune_args = list(),
      enet_tune_args = list(),
      residualizer_args = residualizer_args_none
    ),
    
    # ---- XGB residualisation ---------------------------------------------------------------------
    make_model_spec(
      name = "clpm_xgb_residualized",
      residualizer = "xgb",
      sem_model = "clpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = analyst_order,
      residualizer_exclude = exclude_general,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = xgb_tuning,
      enet_tuning = NULL,
      tune_xgb = tune_xgb,
      tune_enet = FALSE,
      xgb_tune_args = xgb_tune_args,
      enet_tune_args = list(),
      residualizer_args = residualizer_args_xgb
    ),
    
    make_model_spec(
      name = "riclpm_xgb_residualized",
      residualizer = "xgb",
      sem_model = "riclpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = analyst_order,
      residualizer_exclude = exclude_general,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = xgb_tuning,
      enet_tuning = NULL,
      tune_xgb = tune_xgb,
      tune_enet = FALSE,
      xgb_tune_args = xgb_tune_args,
      enet_tune_args = list(),
      residualizer_args = residualizer_args_xgb
    ),
    
    make_model_spec(
      name = "dpm_xgb_residualized",
      residualizer = "xgb",
      sem_model = "dpm",
      sem_c_order = 0,
      sem_exclude = NULL,
      residualizer_c_order = analyst_order,
      residualizer_exclude = exclude_general,
      free_loadings = FALSE,
      bootstrap_B = bootstrap_B,
      xgb_tuning = xgb_tuning,
      enet_tuning = NULL,
      tune_xgb = tune_xgb,
      tune_enet = FALSE,
      xgb_tune_args = xgb_tune_args,
      enet_tune_args = list(),
      residualizer_args = residualizer_args_xgb
    )
  )
}


# helper to run one scenario x sample size cell and save it immediately
run_one_condition <- function(
    N,
    scenario_id,
    scenario_label,
    Delta_list,
    analyst_order = 1L,
    exclude_general = NULL
) {
  
  # build all model specifications for this scenario
  model_specs <- make_model_set(
    analyst_order = analyst_order,
    exclude_general = exclude_general
  )
  
  # run the full simulation workflow
  sim_out <- run_simulation_model_set(
    reps = reps,
    N = N,
    T = T_waves,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    model_specs = model_specs,
    burn_in = burn_in,
    bootstrap_seed = 9000 + N + scenario_id,
    cores = cores,
    base_seed = 100000 + 1000 * scenario_id + N
  )
  
  # extract the combined results frame
  results_df <- sim_out$results
  
  # add scenario metadata for later merging / plotting
  results_df$scenario_id <- scenario_id
  results_df$scenario_label <- scenario_label
  results_df$N <- N
  
  # save one data frame per scenario x sample size combination
  file_out <- file.path(
    out_dir,
    sprintf("scenario_%02d_N%05d.rds", scenario_id, N)
  )
  
  saveRDS(results_df, file_out)
  
  message("Saved: ", file_out)
  
  invisible(results_df)
}

# ---- sample sizes --------------------------------------------------------------------------------

sample_sizes <- c(2000, 1000, 300)

# ---- baseline deltas -----------------------------------------------------------------------------

# scenario 1: constant linear confounding effects only
delta_s1 <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = 0,
  include_2way = FALSE,
  include_3way = FALSE
)

# scenarios 2-4: non-linear confounding with both 2-way and 3-way interactions
# rho_int is the fraction of R2_total assigned to the interaction component
delta_s2plus <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = R2_nl / R2_total,
  include_2way = TRUE,
  include_3way = TRUE
)

# ---- scenario definitions ------------------------------------------------------------------------

# note:
# - the helper scripts can generate additional trajectories, but this master runner
#   intentionally runs only the four final thesis scenarios
# - analyst_order describes what every adjustable model gets to observe
# - analyst_order is fixed to 1 throughout this study: analysts only observe main effects
# - in scenario 4 we additionally omit C1 and C2 from every analyst-side adjustment
scenario_defs <- list(
  list(
    id = 1L,
    label = "constant_linear",
    analyst_order = 1L,
    exclude_general = NULL,
    make_Delta = function() {
      generate_Delta_constant(
        Delta_initial = delta_s1$Delta,
        n_waves = T_waves,
        burn_in = burn_in
      )
    }
  ),
  
  list(
    id = 2L,
    label = "constant_nonlinear",
    analyst_order = 1L,
    exclude_general = NULL,
    make_Delta = function() {
      generate_Delta_constant(
        Delta_initial = delta_s2plus$Delta,
        n_waves = T_waves,
        burn_in = burn_in
      )
    }
  ),
  
  list(
    id = 3L,
    label = "stepwise_nonlinear_rank_stable",
    analyst_order = 1L,
    exclude_general = NULL,
    make_Delta = function() {
      generate_Delta_stepwise(
        Delta_initial = delta_s2plus$Delta,
        n_waves = T_waves,
        burn_in = burn_in,
        step_at = step_at,
        R2_old = R2_total,
        R2_new = R2_new
      )
    }
  ),
  
  list(
    id = 4L,
    label = "stepwise_nonlinear_rank_stable_omit_C1_C2",
    analyst_order = 1L,
    exclude_general = omit_vars,
    make_Delta = function() {
      generate_Delta_stepwise(
        Delta_initial = delta_s2plus$Delta,
        n_waves = T_waves,
        burn_in = burn_in,
        step_at = step_at,
        R2_old = R2_total,
        R2_new = R2_new
      )
    }
  )
)

# ---- run the full grid ---------------------------------------------------------------------------

all_results <- list()
counter <- 1L

n_conditions <- length(sample_sizes) * length(scenario_defs)
condition_idx <- 0L

for (N_now in sample_sizes) {
  for (sc in scenario_defs) {
    
    condition_idx <- condition_idx + 1L
    
    message("------------------------------------------------------------")
    message(
      sprintf(
        "Running scenario %d/%d (N=%d, scenario=%d: %s)",
        condition_idx,
        n_conditions,
        N_now,
        sc$id,
        sc$label
      )
    )
    
    # generate the full Delta trajectory for this scenario
    Delta_list_now <- sc$make_Delta()
    
    # run and save this single scenario x sample size cell
    res_df <- run_one_condition(
      N = N_now,
      scenario_id = sc$id,
      scenario_label = sc$label,
      Delta_list = Delta_list_now,
      analyst_order = sc$analyst_order,
      exclude_general = sc$exclude_general
    )
    
    # also keep it in memory so we can optionally save a combined file
    all_results[[counter]] <- res_df
    counter <- counter + 1L
  }
}

# combined file
all_results_df <- dplyr::bind_rows(all_results)
saveRDS(all_results_df, file.path(out_dir, "all_conditions_combined.rds"))
