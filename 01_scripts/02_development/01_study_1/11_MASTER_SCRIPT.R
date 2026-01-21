# this script serves to call all other scripts and produce the resulting plot in the output folder
# ----------------------------------------------------------------------------------------------------------------
#
# Description:
# This script runs a complete end-to-end example of simulation study 1:
# 1) It loads all required packages and helper functions.
# 2) It samples a baseline delta matrix (D1) describing the effects of time-invariant confounders on X and Y at t = 1.
# 3) It turns D1 into two delta trajectories across T waves:
#    - constant: confounder effects stay the same over time
#    - stepwise: confounder effects increase at a specified wave via an R2 jump
# 4) It defines the true lagged-effect matrix A, where:
#    - autoregressive effects are beta (diagonal of A)
#    - cross-lagged effects are gamma (off-diagonals of A)
# 5) It runs the simulation study across multiple replications, fitting all requested models to each simulated dataset.
# 6) It saves the raw results to disk and creates performance plots (relative bias + RMSE over time) for gamma_XY (X -> Y).
# 7) It saves the plot to the output folder.
# ----------------------------------------------------------------------------------------------------------------

# load here package to manage file paths
library(here)

# set seed for reproducibility
set.seed(1234)

# step 0: load the required packages
source(here("01_scripts", "02_development", "01_study_1", "00_packages.R"))

# step 1: sample a delta-matrix
# get the function to sample the delta-matrix
source(here("01_scripts", "02_development", "01_study_1", "01_delta_sampler.R"))

# sample D1 matrix
D1 <- sample_delta_1(
  k = 3,                                               # number of confounders
  R2_total = 0.15,                                     # total confounder R2 at time t = 1
  min_abs = 0.001,                                     # minimum absolute value for each delta
  max_abs = 0.40,                                      # maximum absolute value for each delta
  R2_nonlin = 0                                        # fraction of R2_total allocated to non-linear terms
)

# step 2: make a delta trajectory from the sampled delta-matrix
# get the function to make the delta trajectory
source(here("01_scripts", "02_development", "01_study_1", "02_delta_trajectory.R"))

# make a constant delta trajectory over 5 time points
D_list_constant <- generate_D_constant(
  D1 = D1,                                             # initial delta matrix
  T  = 5                                               # number of time points
)

# make a stepwise delta trajectory over 5 time points
D_list_stepwise <- generate_D_stepwise(
  D1      = D1,                                        # initial delta matrix
  T       = 5,                                         # number of time points
  step_at = 4,                                         # time point at which to step
  old_R2  = 0.15,                                      # old R2 before the step
  new_R2  = 0.40                                       # new R2 after the step
)

# step 3: define the true matrix
# autoregressive (beta) + cross-lagged (gamma) structure
A <- matrix(c(
  0.20, 0,                         # beta_X and gamma_YX
  0.10, 0.20                       # gamma_XY and beta_Y
), nrow = 2, byrow = TRUE)

# step 4: define the confounder covariance matrix
Psi <- diag(3)                     # uncorrelated confounders with var=1

# step 5: call all the functions that the simulation function (09) needs
source(here("01_scripts", "02_development", "01_study_1", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "01_study_1", "04_lavaan_model_string_builder.R"))
source(here("01_scripts", "02_development", "01_study_1", "05_linear_residualiser.R"))
source(here("01_scripts", "02_development", "01_study_1", "06_model_fitters.R"))
source(here("01_scripts", "02_development", "01_study_1", "07_fit_stat_extractors.R"))
source(here("01_scripts", "02_development", "01_study_1", "08_one_replication_wrapper.R"))

# step 6: run the simulation study function
# get the function to run the full simulation study
source(here("01_scripts", "02_development", "01_study_1", "09_simulation_function.R"))

# run the simulation study
# run a small example simulation study
results_sim <- run_simulation_study(
  reps        = 20,                                    # replications
  N           = 5000,                                  # sample size
  T           = 5,                                     # number of time points
  k           = 3,                                     # number of confounders
  scenarios   = c("constant", "stepwise"),             # D scenarios
  D_scenarios = list(                                  # delta trajectories
    constant = D_list_constant,
    stepwise = D_list_stepwise
  ),
  A           = A,                                     # autoregressive (beta) + cross-lagged (gamma) matrix
  Psi         = Psi,                                   # confounder covariance matrix
  rho_extra   = 0.1,                                   # extra correlation among X_t and Y_t
  models_to_run = c(                                   # models to run
    "clpm",
    "riclpm",
    "dpm",
    "adj",
    "lbca"
  ),
  ci_level    = 0.95,                                  # confidence level

  ###########################################################################
  ### CAUTION: do not set above your machine's available cores - 1!!!
  ### use: parallel::detectCores() - 1 to find out your usable cores
  ###########################################################################

  cores       = parallel::detectCores() - 1,           # number of cores for parallel processing
  base_seed   = 1234                                   # base seed for reproducibility
)

# step 7: save the results 
saveRDS(
  results_sim,
  file = here::here("02_data", "01_research_report","study_1_RR_results.rds")
)

# step 8: produce the plot
# get the plotting function
source(here("01_scripts", "02_development", "01_study_1", "10_plotting.R"))

# produce the plot
p <- plot_sim_study_results(
  results_sim = results_sim,
  true_A      = A
)

# print the plot
print(p$combined_gamma_XY)

# save the plot
ggsave(
  filename = here::here("03_output", "01_latest", "RR_study_1.png"),
  plot     = p$combined_gamma_XY,
  width    = 10,
  height   = 8,
  dpi      = 300
)
