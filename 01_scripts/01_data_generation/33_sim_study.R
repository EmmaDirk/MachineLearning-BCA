# =========================================================================================
# sim_study_reduced_s1235.R
# 
# Reduced simulation study.
# 
# We always fit the requested reduced model set to the generated scenario.
# 
# We always use:
# - T = 5
# - burn_in = 20
# - R^2_total = 0.05 
# - k = 5
# - Phi = {0.25, 0.00
#          0.10, 0.25}
# - Sigma = {1.00, 0.30
#            0.30, 1.00}
# 
# We test the following scenarios:
# 
# A) 2000 observations
# 1) constant linear confounding effects
# 2) constant non-linear confounding effects, R^2_nl = 0.02
# 3) stepwise rank-stable non-linear confounding effects,
#    R^2_nl = 0.02, new_R^2 = 0.10, step_at = 3
# 5) the same specification as 3, but we omit c1 and c2 
#    from every analyst-side adjustment.
#
# B) 1000 observations: scenarios 1, 2, 3, 5
#
# C) 300 observations: scenarios 1, 2, 3, 5
#
# For every scenario we fit:
# - the CLPM model with no adjustment
# - the RI-CLPM model with no adjustment
# - the DPM model with no adjustment
#
# - the CLPM model with observed linear confounder control
# - the RI-CLPM model with observed linear confounder control
# - the DPM model with observed linear confounder control
#
# - the CLPM model after XGB residualisation
# - the RI-CLPM model after XGB residualisation
# - the DPM model after XGB residualisation
#
# In the non-linear scenarios we use both 2-way and 3-way interactions in the DGM.
# The analyst-side linear confounder models use only main effects.
#
# We save one data frame per scenario x sample size combination.
#
# =========================================================================================

# --------------------------------------- logistics ----------------------------------------

library(here)

# reusable directories
script_dir <- here::here("Simulation Studies", "01_thesis_study_v1", "01_scripts")
out_dir    <- here::here("Simulation Studies", "01_thesis_study_v1", "03_output")

# load the study scripts
# Keep one explicit master list and one worker list so the cluster uses the exact
# same code versions as the master session during the replication stage.
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

for (f in simulation_script_files) {
  source(file.path(script_dir, f))
}

# reproducibility
set.seed(1233)

# IMPORTANT:
# this study parallelises across replications at the R-process level.
# Therefore we force threaded BLAS / OpenMP libraries to use one thread per R process.
# Otherwise a large server can become slower than a laptop due to thread oversubscription.
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

# --------------------------------------- shared DGM ----------------------------------------

# number of confounders
k <- 5

# number of waves
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
lambda_mix <- 0.5

# analyst-side omissions used only in scenario 5
omit_vars <- c("c1", "c2")

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

# --------------------------------------- shared model specs --------------------------------

# global 
reps <- 2000
bootstrap_B <- 40

# IMPORTANT:
# do not start with all available server cores.
# After capping BLAS / OpenMP threads to 1, a moderate worker count is usually faster and
# much more stable than 100+ PSOCK workers for this workload. Start here, then scale up.
cores <- 50

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
  nthread = 32,
  seed = 123
)

residualizer_args_xgb <- list(
  oof_folds = 8,
  use_oof = TRUE,
  nthread = 1,
  seed = 123
)

# elastic net 
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

# --------------------------------------- helper functions ----------------------------------

# helper to build the reduced model set for one scenario
# analyst_order describes what every adjustable model gets to observe
# exclude_general applies to every adjustable model
make_model_set <- function(
    analyst_order = 1L,
    exclude_general = NULL
) {
  
  list(
    # ----------------------------- no adjustment models -----------------------------
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
    
    # ----------------------- observed linear confounder control ----------------------
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
    
    # ----------------------------- XGB residualisation ------------------------------
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
    n_waves = T_waves,
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
    sprintf("reduced_s%02d_N%05d.rds", scenario_id, N)
  )
  
  saveRDS(results_df, file_out)
  
  message("Saved: ", file_out)
  
  invisible(results_df)
}

# --------------------------------------- sample sizes --------------------------------------

sample_sizes <- c(2000, 1000, 300)

# --------------------------------------- baseline deltas ----------------------------------

# scenario 1: constant linear confounding effects only
delta_s1 <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = 0,
  include_2way = FALSE,
  include_3way = FALSE
)

# scenarios 2-6: non-linear confounding with both 2-way and 3-way interactions
# rho_int is the fraction of total R2 assigned to the nonlinear component
delta_s2plus <- sample_delta_t(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  rho_int = R2_nl / R2_total,
  include_2way = TRUE,
  include_3way = TRUE
)

# --------------------------------------- scenario definitions ------------------------------

# note:
# - analyst_order describes what every adjustable model gets to observe
# - analyst_order is fixed to 1 throughout this study: analysts only observe main effects
# - in scenario 5 we additionally omit c1 and c2 from every analyst-side adjustment
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
    id = 5L,
    label = "stepwise_nonlinear_rank_stable_omit_c1_c2",
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

# --------------------------------------- run the full grid --------------------------------

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

# convenience combined file
all_results_df <- dplyr::bind_rows(all_results)
saveRDS(all_results_df, file.path(out_dir, "reduced_all_conditions_combined.rds"))
