# This script is a compact end-to-end example of the new engine.
#
# The example is deliberately small.
# Its purpose is not to be the final study.
# Its purpose is to show the intended workflow:
#
# 1) source all scripts
# 2) sample the baseline delta matrix
# 3) turn it into one or more delta trajectories
# 4) define Phi, Omega11, and Sigma
# 5) fit the model functions directly on one simulated data set
# 6) optionally build a model_specs list and run the simulation study
# ------------------------------------------------------------------------------------------

# source everything
source("11_source_all.R")

# optional: check the packages you will need
# if you want to use the xgb residualiser, set include_xgboost = TRUE
check_sim_engine_packages(include_xgboost = FALSE)

# ------------------------- design choices -------------------------

set.seed(1234)

# number of confounders
k <- 3

# number of waves
T <- 5

# covariance matrix of the base confounders
Omega11 <- diag(k)

# target covariance matrix of (X_t, Y_t)
Sigma <- matrix(
  c(1.00, 0.30,
    0.30, 1.00),
  nrow = 2,
  byrow = TRUE
)

# lag matrix:
# diagonal   = autoregressive effects
# off-diagonals = cross-lagged effects
Phi <- matrix(
  c(0.20, 0.00,
    0.10, 0.20),
  nrow = 2,
  byrow = TRUE
)

# ------------------------- delta setup -------------------------

# sample the baseline confounder-effect matrix at wave 1
Delta1 <- sample_delta_1(
  k = k,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0.00,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0.001,
  max_abs = 0.40
)

# turn that baseline matrix into two trajectories
Delta_constant <- generate_Delta_constant(
  Delta1 = Delta1,
  T = T
)

Delta_stepwise <- generate_Delta_stepwise(
  Delta1 = Delta1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

Delta_scenarios <- list(
  constant = Delta_constant,
  stepwise = Delta_stepwise
)

# ------------------------- one data set -------------------------

df <- simulate_panel_data(
  N = 300,
  T = T,
  Phi = Phi,
  Delta_list = Delta_constant,
  Omega11 = Omega11,
  Sigma = Sigma,
  return_full = FALSE,
  seed = 999
)

# ------------------------- direct model calls -------------------------
# this is the main modelling workflow now:
# you can just call the fitting functions directly.

fit1 <- fit_clpm(
  df = df,
  T = T,
  residualiser = "none"
)

fit2 <- fit_riclpm(
  df = df,
  T = T,
  residualiser = "linear",
  resid_k = k,
  resid_interaction_order = 1,
  bootstrap_R = 99
)

fit3 <- fit_dpm(
  df = df,
  T = T,
  residualiser = "none",
  free_loadings = TRUE
)

# inspect the minimal returned objects
fit1$converged
fit1$proper
fit1$bic
fit1$parameters

# ------------------------- full simulation example -------------------------
# model_specs tells the runner which models to call in each replication.

model_specs <- list(
  clpm = list(
    fun = fit_clpm,
    args = list(
      T = T,
      residualiser = "none"
    )
  ),

  riclpm_linear = list(
    fun = fit_riclpm,
    args = list(
      T = T,
      residualiser = "linear",
      resid_k = k,
      resid_interaction_order = 1,
      bootstrap_R = 49
    )
  ),

  dpm_free = list(
    fun = fit_dpm,
    args = list(
      T = T,
      residualiser = "none",
      free_loadings = TRUE
    )
  )
)

# small example run
results_sim <- run_simulation_study(
  reps = 20,
  N = 300,
  T = T,
  Delta_scenarios = Delta_scenarios,
  Phi = Phi,
  Omega11 = Omega11,
  Sigma = Sigma,
  model_specs = model_specs,
  cores = 1,
  base_seed = 1234
)

head(results_sim)
