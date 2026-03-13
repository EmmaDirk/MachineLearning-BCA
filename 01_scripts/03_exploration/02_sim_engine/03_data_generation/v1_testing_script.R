# libraries
library(here)
library(ggplot2)
library(viridis)

# source scripts
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "02_delta_trajectory", "v1_delta_trajectory.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "03_data_generation", "v1_data_generation.R"))

# reproducibility
set.seed(123)

# create an omega_11
Omega11 <- matrix(c(
  1, 0.3, 0.2,
  0.3, 1, 0.4,
  0.2, 0.4, 1
), nrow = 3, byrow = TRUE)

# sample delta
out1 <- sample_delta_1(
  k = 3,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0.3,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

# extract delta
delta1 <- out1$Delta
Omega <- out1$Omega

# create three scenarios:
# constant
Delta_const <- generate_Delta_constant(
  Delta1 = delta1,
  T  = 5
)

# stepwise
Delta_step <- generate_Delta_stepwise(
  Delta1 = delta1,
  T = 5,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# stepwise mixture
Delta_mix <- generate_Delta_stepwise_mixture(
  Delta1 = delta1,
  T = 5,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.3,
  lambda_NL = 0.3,

  # important: this seed must differ from the seed used above for sampling Delta1
  seed = 1234
)

# look at the results
# visually inspect
Delta_const
Delta_step
Delta_mix

# phi matrix (2x2)
phi <- matrix(c(
  0.3, 0,
  0.1, 0.3
), nrow = 2, byrow = TRUE)

# sigma matrix (2x2)
Sigma <- matrix(c(
  1, 0.3,
  0.3, 1
), nrow = 2, byrow = TRUE)

# generate data using the constant scenario
out_const <- simulate_panel_data(
  N = 1000000,
  T = 5,
  Phi = phi,
  Omega11 = Omega11,
  Sigma = Sigma,
  Delta_list = Delta_const,
  eig_tol = 1e-10
)

data_const <- out_const$data
cov(data_const)
