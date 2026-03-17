# to test the fit functions
# this is to later be expanded to also deep dive into all the helpers and wrappers
# we take the following steps:
# 1) generate some testing data
# 2) run the CLPM fitter with various options
# 3) run the DPM fitter with various options
# 4) run the RI-CLPM fitter with various options

set.seed(1234)

# ------------------------------- 1) generate testing data -------------------------------
# step 1.1: load the required packages
library(here)
library(tidyverse)
library(mvtnorm)
library(xgboost)
library(knitr)

# step 1.2: load the required functions for data simulation
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "02_delta_trajectory", "v1_delta_trajectory.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "03_data_generation", "v1_data_generation.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "04_models", "v1_residualisers.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "04_models", "v1_lavaan.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "05_logistics", "v1_fitter.R"))

# step 1.3: define some parameters
T <- 5
N <- 5000

Sigma <- matrix(c(
  1.0, 0.30,
  0.30, 1.0
), nrow = 2, byrow = TRUE)

Phi <- matrix(c(
  0.40, 0,
  0.10, 0.35
), nrow = 2, byrow = TRUE)

Omega11 <- matrix(c(
  1.0, 0.3, 0.2,
  0.3, 1.0, 0.4,
  0.2, 0.4, 1.0
), nrow = 3, byrow = TRUE)

# step 1.4: sample delta coefficients
out_delta <- sample_delta_1(
  k = 3,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0.1,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

Delta <- out_delta$Delta

# step 1.5: generate a constant scenario
delta_list <- generate_Delta_constant(
  Delta1 = Delta,
  T = T
)

# step 1.6: generate the data
out <- simulate_panel_data(
  N = N,
  T = T,
  Phi = Phi,
  Delta_list = delta_list,
  Omega11 = Omega11,
  Sigma = Sigma
)

data <- out$data 

# -------------------------------- 2) run the CLPM fitter with various options --------------------------------

# 2.1: run the CLPM fitter with no residualiser
fit_clpm_none <- fit_clpm(
  df = data,
  T = T,
  residualiser = "none",

  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  model_k = 0,
  model_exclude = NULL,
  model_confounder_order = 0
)

# 2.2: run the CLPM fitter with linear residualiser
fit_clpm_linear <- fit_clpm(
  df = data,
  T = T,
  residualiser = "linear",

  resid_k = 1,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  model_k = 0,
  model_exclude = NULL,
  model_confounder_order = 0,

  bootstrap_R = 200,
  bootstrap_seed = 123
)
