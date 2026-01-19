# this script serves to call all other scripts and produce the resulting plot in the output folder
# ----------------------------------------------------------------------------------------------------------------

# load here package to manage file paths
library(here)

# set seed for reproducibility
set.seed(1234)

# step 0: load the required packages
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "00_packages.R"))

# step 1: sample a beta-matrix
# get the function to sample the beta-matrix
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "01_beta_sampler.R"))

# sample B1 matrix
B1 <- sample_B1(
  k = 3,                                               # number of confounders
  R2_1 = 0.15,                                         # total confounder R2 at time t = 1 
  min_abs = 0.001,                                     # minimum absolute value for each beta
  max_abs = 0.40,                                      # maximum absolute value for each beta
  eta_1 = 0                                            # fraction of R2_1 allocated to non-linear terms
)

# step 2: make a beta trajectory from the sampled beta-matrix
# get the function to make the beta trajectory
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "02_beta_trajectory.R"))

# make a constant beta trajectory over 5 time points
B_list_constant <- generate_B_constant(
  B1 = B1,                                             # initial beta matrix
  T = 5                                                # number of time points
)

# make a stepwise beta trajectory over 5 time points
B_list_stepwise <- generate_B_stepwise(
  B1     = B1,                                         # initial beta matrix
  T      = 5,                                          # number of time points
  step_at = 4,                                         # time point at which to step
  old_R2 = 0.15,                                       # old R2 before the step
  new_R2 = 0.40                                        # new R2 after the step
)

# step 3: define the true matrix
# autoregressive + cross-lag structure
A <- matrix(c(
  0.20, 0,                         # autoregressive and cross-lag for Y
  0.10, 0.20                       # cross-lag and autoregressive for X
), nrow = 2, byrow = TRUE)

# step 4: define the confounder covariance matrix
Psi <- diag(3)                     # uncorrelated confounders with var=1

# step 5: call all the functions that the simulation function (09) needs
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "04_lavaan_model_string_builder.R"))
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "05_linear_residualiser.R"))
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "06_model_fitters.R"))
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "07_fit_stat_extractors.R"))
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "08_one_replication_wrapper.R"))

# step 6: run the simulation study function
# get the function to run the full simulation study
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "09_simulation_function.R"))

# run the simulation study
# run a small example simulation study
results_sim <- run_simulation_study(
  reps        = 20,                                    # replications
  N           = 5000,                                  # sample size
  T           = 5,                                     # number of time points
  k           = 3,                                     # number of confounders
  scenarios   = c("constant", "stepwise"),             # B scenarios
  B_scenarios = list(                                  # beta trajectories
    constant = B_list_constant,
    stepwise = B_list_stepwise
  ),
  A           = A,                                     # autoregressive + cross-lag matrix
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
source(here("01_scripts", "02_development", "00_research_report_engines", "01_study_1", "10_plotting.R"))

# produce the plot
p <- plot_sim_study_results(
  results_sim = results_sim,
  true_A      = A
)

# print the plot
print(p$combined_XY)

# save the plot
ggsave(
  filename = here::here("03_output", "01_latest", "RR_study_1.png"),
  plot     = p$combined_XY,
  width    = 10,
  height   = 8,
  dpi      = 300
)
