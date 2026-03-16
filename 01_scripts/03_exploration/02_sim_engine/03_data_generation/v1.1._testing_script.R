# goal is to stress test all parts of the panel-data simulation process
# to make sure that everything is working as expected.
# focus:
# - covariance of confounder features
# - correctness of Psi_t construction
# - correctness of M_t recursion
# - correctness of the state variance equation in vacuum
# - correctness of implied covariance of (X_t, Y_t)
# - behavior under k = 1, 2, 3, 4 and different Delta trajectories

###########################################################################
# PANEL DATA SIMULATION
###########################################################################

library(here)
library(mvtnorm)

source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "02_delta_trajectory", "v1_delta_trajectory.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "03_panel_sim", "v1_simulate_panel_data.R"))

###########################################################################
# GLOBAL SETUP
###########################################################################

N_big <- 200000
N_small <- 5000
T_panel <- 5

Sigma <- matrix(c(
  1.0, 0.30,
  0.30, 1.0
), nrow = 2, byrow = TRUE)

Phi <- matrix(c(
  0.40, 0.15,
  0.10, 0.35
), nrow = 2, byrow = TRUE)

Sigma
Phi

###########################################################################
# CASE 1: k = 1, NO INTERACTIONS
###########################################################################

Omega11_1k <- diag(1)

base_1 <- sample_delta_1(
  k = 1,
  Omega11 = Omega11_1k,
  R2_total = 0.15,
  R2_interaction = 0,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

Delta1_1 <- base_1$Delta
Delta_list_1_const <- generate_Delta_constant(
  Delta1 = Delta1_1,
  T = T_panel
)

Delta1_1
Delta_list_1_const

sim_1 <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_1_const,
  Omega11 = Omega11_1k,
  Sigma = Sigma,
  seed = 123
)

names(sim_1)

# check base confounders
head(sim_1$base_confounders)
cov(sim_1$base_confounders)

# for k = 1 the full confounder matrix should just be the base confounder
head(sim_1$confounders)
cov(sim_1$confounders)
sim_1$Omega_full

# compare empirical confounder covariance to analytic covariance
cov(sim_1$confounders)
sim_1$Omega_full
cov(sim_1$confounders) - sim_1$Omega_full

# check wave-specific covariance of observed variables
cov(sim_1$data[, c("x1", "y1")])
cov(sim_1$data[, c("x2", "y2")])
cov(sim_1$data[, c("x3", "y3")])
cov(sim_1$data[, c("x4", "y4")])
cov(sim_1$data[, c("x5", "y5")])

# compare to target Sigma
cov(sim_1$data[, c("x1", "y1")]) - Sigma
cov(sim_1$data[, c("x2", "y2")]) - Sigma
cov(sim_1$data[, c("x3", "y3")]) - Sigma
cov(sim_1$data[, c("x4", "y4")]) - Sigma
cov(sim_1$data[, c("x5", "y5")]) - Sigma

# inspect Psi matrices
sim_1$Psi_list[[1]]
sim_1$Psi_list[[2]]
sim_1$Psi_list[[3]]
sim_1$Psi_list[[4]]
sim_1$Psi_list[[5]]

# inspect eigenvalues of Psi matrices
eigen(sim_1$Psi_list[[1]], symmetric = TRUE)$values
eigen(sim_1$Psi_list[[2]], symmetric = TRUE)$values
eigen(sim_1$Psi_list[[3]], symmetric = TRUE)$values
eigen(sim_1$Psi_list[[4]], symmetric = TRUE)$values
eigen(sim_1$Psi_list[[5]], symmetric = TRUE)$values

# inspect M recursion objects
sim_1$M_list[[1]]
sim_1$M_list[[2]]
sim_1$M_list[[3]]
sim_1$M_list[[4]]
sim_1$M_list[[5]]

###########################################################################
# CASE 2: k = 2, TWO-WAY INTERACTION ONLY
###########################################################################

Omega11_2 <- diag(2)

base_2 <- sample_delta_1(
  k = 2,
  Omega11 = Omega11_2,
  R2_total = 0.15,
  R2_interaction = 0.30,
  include_2way = TRUE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

Delta1_2 <- base_2$Delta
Delta_list_2_const <- generate_Delta_constant(
  Delta1 = Delta1_2,
  T = T_panel
)

Delta1_2
colnames(Delta1_2)

sim_2 <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_2_const,
  Omega11 = Omega11_2,
  Sigma = Sigma,
  seed = 123
)

# check feature covariance analytically vs empirically
colnames(sim_2$confounders)
head(sim_2$confounders)

cov(sim_2$confounders)
sim_2$Omega_full
cov(sim_2$confounders) - sim_2$Omega_full

# for k = 2 with independent base confounders, expected structure:
# Cov(c1, c2) = 0
# Cov(c1, c1:2) = 0
# Cov(c2, c1:2) = 0
# variances all 1
sim_2$Omega_full

# observed covariance of (X_t, Y_t) at each wave
cov(sim_2$data[, c("x1", "y1")])
cov(sim_2$data[, c("x2", "y2")])
cov(sim_2$data[, c("x3", "y3")])
cov(sim_2$data[, c("x4", "y4")])
cov(sim_2$data[, c("x5", "y5")])

cov(sim_2$data[, c("x1", "y1")]) - Sigma
cov(sim_2$data[, c("x2", "y2")]) - Sigma
cov(sim_2$data[, c("x3", "y3")]) - Sigma
cov(sim_2$data[, c("x4", "y4")]) - Sigma
cov(sim_2$data[, c("x5", "y5")]) - Sigma

# inspect Psi and M
sim_2$Psi_list[[1]]
sim_2$Psi_list[[2]]
sim_2$Psi_list[[3]]

eigen(sim_2$Psi_list[[1]], symmetric = TRUE)$values
eigen(sim_2$Psi_list[[2]], symmetric = TRUE)$values
eigen(sim_2$Psi_list[[3]], symmetric = TRUE)$values

sim_2$M_list[[1]]
sim_2$M_list[[2]]
sim_2$M_list[[3]]

###########################################################################
# CASE 3: k = 3, INDEPENDENT BASE CONFOUNDERS, TWO-WAY + THREE-WAY
###########################################################################

Omega11_3 <- diag(3)

base_3 <- sample_delta_1(
  k = 3,
  Omega11 = Omega11_3,
  R2_total = 0.15,
  R2_interaction = 0.30,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

Delta1_3 <- base_3$Delta
Delta_list_3_const <- generate_Delta_constant(
  Delta1 = Delta1_3,
  T = T_panel
)

Delta1_3
colnames(Delta1_3)

sim_3 <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_3_const,
  Omega11 = Omega11_3,
  Sigma = Sigma,
  seed = 123
)

# inspect analytic full Omega for features
sim_3$Omega_full

# inspect empirical confounder covariance
cov(sim_3$confounders)

# compare empirical to analytic
cov(sim_3$confounders) - sim_3$Omega_full

# observed covariance by wave
cov(sim_3$data[, c("x1", "y1")])
cov(sim_3$data[, c("x2", "y2")])
cov(sim_3$data[, c("x3", "y3")])
cov(sim_3$data[, c("x4", "y4")])
cov(sim_3$data[, c("x5", "y5")])

cov(sim_3$data[, c("x1", "y1")]) - Sigma
cov(sim_3$data[, c("x2", "y2")]) - Sigma
cov(sim_3$data[, c("x3", "y3")]) - Sigma
cov(sim_3$data[, c("x4", "y4")]) - Sigma
cov(sim_3$data[, c("x5", "y5")]) - Sigma

# inspect Psi and eigenvalues
sim_3$Psi_list[[1]]
sim_3$Psi_list[[2]]
sim_3$Psi_list[[3]]
sim_3$Psi_list[[4]]
sim_3$Psi_list[[5]]

eigen(sim_3$Psi_list[[1]], symmetric = TRUE)$values
eigen(sim_3$Psi_list[[2]], symmetric = TRUE)$values
eigen(sim_3$Psi_list[[3]], symmetric = TRUE)$values
eigen(sim_3$Psi_list[[4]], symmetric = TRUE)$values
eigen(sim_3$Psi_list[[5]], symmetric = TRUE)$values

# inspect M recursion
sim_3$M_list[[1]]
sim_3$M_list[[2]]
sim_3$M_list[[3]]
sim_3$M_list[[4]]
sim_3$M_list[[5]]

###########################################################################
# CASE 4: k = 4, CORRELATED BASE CONFOUNDERS, TWO-WAY + THREE-WAY
###########################################################################

Omega11_4 <- matrix(c(
  1.0, 0.3, 0.2, 0.1,
  0.3, 1.0, 0.4, 0.2,
  0.2, 0.4, 1.0, 0.3,
  0.1, 0.2, 0.3, 1.0
), nrow = 4, byrow = TRUE)

base_4 <- sample_delta_1(
  k = 4,
  Omega11 = Omega11_4,
  R2_total = 0.15,
  R2_interaction = 0.30,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

Delta1_4 <- base_4$Delta
Delta_list_4_const <- generate_Delta_constant(
  Delta1 = Delta1_4,
  T = T_panel
)

Delta1_4
colnames(Delta1_4)

sim_4 <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 123
)

# inspect analytic Omega_full
sim_4$Omega_full

# empirical feature covariance
cov(sim_4$confounders)

# difference empirical - analytic
cov(sim_4$confounders) - sim_4$Omega_full

# wave-specific covariance of (X_t, Y_t)
cov(sim_4$data[, c("x1", "y1")])
cov(sim_4$data[, c("x2", "y2")])
cov(sim_4$data[, c("x3", "y3")])
cov(sim_4$data[, c("x4", "y4")])
cov(sim_4$data[, c("x5", "y5")])

cov(sim_4$data[, c("x1", "y1")]) - Sigma
cov(sim_4$data[, c("x2", "y2")]) - Sigma
cov(sim_4$data[, c("x3", "y3")]) - Sigma
cov(sim_4$data[, c("x4", "y4")]) - Sigma
cov(sim_4$data[, c("x5", "y5")]) - Sigma

# inspect Psi matrices
sim_4$Psi_list[[1]]
sim_4$Psi_list[[2]]
sim_4$Psi_list[[3]]
sim_4$Psi_list[[4]]
sim_4$Psi_list[[5]]

# check positive semidefiniteness
eigen(sim_4$Psi_list[[1]], symmetric = TRUE)$values
eigen(sim_4$Psi_list[[2]], symmetric = TRUE)$values
eigen(sim_4$Psi_list[[3]], symmetric = TRUE)$values
eigen(sim_4$Psi_list[[4]], symmetric = TRUE)$values
eigen(sim_4$Psi_list[[5]], symmetric = TRUE)$values

# inspect M recursion
sim_4$M_list[[1]]
sim_4$M_list[[2]]
sim_4$M_list[[3]]
sim_4$M_list[[4]]
sim_4$M_list[[5]]

###########################################################################
# CHECK Psi_1 AND M_1 IN A VACUUM FOR k = 4
###########################################################################

Delta1_test <- Delta_list_4_const[[1]]
Omega_full_test <- sim_4$Omega_full

Delta1_test
Omega_full_test

# theoretical M1
M1_manual <- Delta1_test %*% Omega_full_test
M1_manual

# compare to returned M1
sim_4$M_list[[1]]
M1_manual - sim_4$M_list[[1]]

# theoretical Psi1
Psi1_manual <- Sigma - Delta1_test %*% Omega_full_test %*% t(Delta1_test)
Psi1_manual

# compare to returned Psi1
sim_4$Psi_list[[1]]
Psi1_manual - sim_4$Psi_list[[1]]

###########################################################################
# CHECK Psi_t AND M_t RECURSION IN A VACUUM FOR k = 4
###########################################################################

Delta2_test <- Delta_list_4_const[[2]]
M1_test <- sim_4$M_list[[1]]

# manual M2
M2_manual <- Phi %*% M1_test + Delta2_test %*% Omega_full_test
M2_manual

# returned M2
sim_4$M_list[[2]]
M2_manual - sim_4$M_list[[2]]

# manual Psi2
Psi2_manual <-
  Sigma -
  Phi %*% Sigma %*% t(Phi) -
  Delta2_test %*% Omega_full_test %*% t(Delta2_test) -
  Phi %*% M1_test %*% t(Delta2_test) -
  Delta2_test %*% t(M1_test) %*% t(Phi)

Psi2_manual

# returned Psi2
sim_4$Psi_list[[2]]
Psi2_manual - sim_4$Psi_list[[2]]

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: WAVE 1
###########################################################################

Delta1_test <- Delta_list_4_const[[1]]
Omega_full_test <- sim_4$Omega_full
Psi1_test <- sim_4$Psi_list[[1]]

Delta1_test
Omega_full_test
Psi1_test

# reconstruct Var(W1) directly from the model equation
Var_W1_manual <- Delta1_test %*% Omega_full_test %*% t(Delta1_test) + Psi1_test
Var_W1_manual

# compare to target Sigma
Sigma
Var_W1_manual - Sigma

# same thing, written exactly as in the derivation
Delta1_test %*% Omega_full_test %*% t(Delta1_test)
Psi1_test
Delta1_test %*% Omega_full_test %*% t(Delta1_test) + Psi1_test

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: WAVE 2
###########################################################################

Delta2_test <- Delta_list_4_const[[2]]
M1_test <- sim_4$M_list[[1]]
Psi2_test <- sim_4$Psi_list[[2]]

Delta2_test
M1_test
Psi2_test

# reconstruct Var(W2) term by term
term_AR_2 <- Phi %*% Sigma %*% t(Phi)
term_D_2 <- Delta2_test %*% Omega_full_test %*% t(Delta2_test)
term_cross_2_a <- Phi %*% M1_test %*% t(Delta2_test)
term_cross_2_b <- Delta2_test %*% t(M1_test) %*% t(Phi)

term_AR_2
term_D_2
term_cross_2_a
term_cross_2_b
Psi2_test

Var_W2_manual <- term_AR_2 + term_D_2 + term_cross_2_a + term_cross_2_b + Psi2_test
Var_W2_manual

# compare to target Sigma
Sigma
Var_W2_manual - Sigma

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: WAVES 3, 4, 5
###########################################################################

# wave 3
Delta3_test <- Delta_list_4_const[[3]]
M2_test <- sim_4$M_list[[2]]
Psi3_test <- sim_4$Psi_list[[3]]

Var_W3_manual <-
  Phi %*% Sigma %*% t(Phi) +
  Delta3_test %*% Omega_full_test %*% t(Delta3_test) +
  Phi %*% M2_test %*% t(Delta3_test) +
  Delta3_test %*% t(M2_test) %*% t(Phi) +
  Psi3_test

Var_W3_manual
Var_W3_manual - Sigma

# wave 4
Delta4_test <- Delta_list_4_const[[4]]
M3_test <- sim_4$M_list[[3]]
Psi4_test <- sim_4$Psi_list[[4]]

Var_W4_manual <-
  Phi %*% Sigma %*% t(Phi) +
  Delta4_test %*% Omega_full_test %*% t(Delta4_test) +
  Phi %*% M3_test %*% t(Delta4_test) +
  Delta4_test %*% t(M3_test) %*% t(Phi) +
  Psi4_test

Var_W4_manual
Var_W4_manual - Sigma

# wave 5
Delta5_test <- Delta_list_4_const[[5]]
M4_test <- sim_4$M_list[[4]]
Psi5_test <- sim_4$Psi_list[[5]]

Var_W5_manual <-
  Phi %*% Sigma %*% t(Phi) +
  Delta5_test %*% Omega_full_test %*% t(Delta5_test) +
  Phi %*% M4_test %*% t(Delta5_test) +
  Delta5_test %*% t(M4_test) %*% t(Phi) +
  Psi5_test

Var_W5_manual
Var_W5_manual - Sigma

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: ALL WAVES IN ONE LOOP
###########################################################################

for (t in 1:T_panel) {

  cat("\n==============================\n")
  cat("wave =", t, "\n")
  cat("==============================\n")

  Delta_t <- Delta_list_4_const[[t]]
  Psi_t <- sim_4$Psi_list[[t]]

  if (t == 1) {

    Var_W_t <- Delta_t %*% sim_4$Omega_full %*% t(Delta_t) + Psi_t

  } else {

    M_prev <- sim_4$M_list[[t - 1]]

    Var_W_t <-
      Phi %*% Sigma %*% t(Phi) +
      Delta_t %*% sim_4$Omega_full %*% t(Delta_t) +
      Phi %*% M_prev %*% t(Delta_t) +
      Delta_t %*% t(M_prev) %*% t(Phi) +
      Psi_t
  }

  print(Var_W_t)
  print(Var_W_t - Sigma)
}

###########################################################################
# STEPWISE Delta TRAJECTORY
###########################################################################

Delta_list_4_step <- generate_Delta_stepwise(
  Delta1 = Delta1_4,
  T = T_panel,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.30
)

Delta_list_4_step

sim_4_step <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_step,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 123
)

# observed covariance by wave should still be Sigma
cov(sim_4_step$data[, c("x1", "y1")])
cov(sim_4_step$data[, c("x2", "y2")])
cov(sim_4_step$data[, c("x3", "y3")])
cov(sim_4_step$data[, c("x4", "y4")])
cov(sim_4_step$data[, c("x5", "y5")])

cov(sim_4_step$data[, c("x1", "y1")]) - Sigma
cov(sim_4_step$data[, c("x2", "y2")]) - Sigma
cov(sim_4_step$data[, c("x3", "y3")]) - Sigma
cov(sim_4_step$data[, c("x4", "y4")]) - Sigma
cov(sim_4_step$data[, c("x5", "y5")]) - Sigma

# inspect change in Psi around the step
sim_4_step$Psi_list[[1]]
sim_4_step$Psi_list[[2]]
sim_4_step$Psi_list[[3]]
sim_4_step$Psi_list[[4]]
sim_4_step$Psi_list[[5]]

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: STEPWISE TRAJECTORY
###########################################################################

Omega_full_step <- sim_4_step$Omega_full

# wave 1
Delta1_step <- Delta_list_4_step[[1]]
Psi1_step <- sim_4_step$Psi_list[[1]]

Var_W1_step <- Delta1_step %*% Omega_full_step %*% t(Delta1_step) + Psi1_step
Var_W1_step
Var_W1_step - Sigma

# wave 2
Delta2_step <- Delta_list_4_step[[2]]
M1_step <- sim_4_step$M_list[[1]]
Psi2_step <- sim_4_step$Psi_list[[2]]

Var_W2_step <-
  Phi %*% Sigma %*% t(Phi) +
  Delta2_step %*% Omega_full_step %*% t(Delta2_step) +
  Phi %*% M1_step %*% t(Delta2_step) +
  Delta2_step %*% t(M1_step) %*% t(Phi) +
  Psi2_step

Var_W2_step
Var_W2_step - Sigma

# wave 3
Delta3_step <- Delta_list_4_step[[3]]
M2_step <- sim_4_step$M_list[[2]]
Psi3_step <- sim_4_step$Psi_list[[3]]

Var_W3_step <-
  Phi %*% Sigma %*% t(Phi) +
  Delta3_step %*% Omega_full_step %*% t(Delta3_step) +
  Phi %*% M2_step %*% t(Delta3_step) +
  Delta3_step %*% t(M2_step) %*% t(Phi) +
  Psi3_step

Var_W3_step
Var_W3_step - Sigma

# wave 4
Delta4_step <- Delta_list_4_step[[4]]
M3_step <- sim_4_step$M_list[[3]]
Psi4_step <- sim_4_step$Psi_list[[4]]

Var_W4_step <-
  Phi %*% Sigma %*% t(Phi) +
  Delta4_step %*% Omega_full_step %*% t(Delta4_step) +
  Phi %*% M3_step %*% t(Delta4_step) +
  Delta4_step %*% t(M3_step) %*% t(Phi) +
  Psi4_step

Var_W4_step
Var_W4_step - Sigma

# wave 5
Delta5_step <- Delta_list_4_step[[5]]
M4_step <- sim_4_step$M_list[[4]]
Psi5_step <- sim_4_step$Psi_list[[5]]

Var_W5_step <-
  Phi %*% Sigma %*% t(Phi) +
  Delta5_step %*% Omega_full_step %*% t(Delta5_step) +
  Phi %*% M4_step %*% t(Delta5_step) +
  Delta5_step %*% t(M4_step) %*% t(Phi) +
  Psi5_step

Var_W5_step
Var_W5_step - Sigma

###########################################################################
# MIXTURE Delta TRAJECTORY
###########################################################################

Delta_list_4_mix <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_4,
  T = T_panel,
  Omega = base_4$Omega,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.30,
  lambda_L = 0.50,
  lambda_NL = 0.50,
  seed = 123
)

Delta_list_4_mix

sim_4_mix <- simulate_panel_data(
  N = N_big,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_mix,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 123
)

# observed covariance by wave should still be Sigma
cov(sim_4_mix$data[, c("x1", "y1")])
cov(sim_4_mix$data[, c("x2", "y2")])
cov(sim_4_mix$data[, c("x3", "y3")])
cov(sim_4_mix$data[, c("x4", "y4")])
cov(sim_4_mix$data[, c("x5", "y5")])

cov(sim_4_mix$data[, c("x1", "y1")]) - Sigma
cov(sim_4_mix$data[, c("x2", "y2")]) - Sigma
cov(sim_4_mix$data[, c("x3", "y3")]) - Sigma
cov(sim_4_mix$data[, c("x4", "y4")]) - Sigma
cov(sim_4_mix$data[, c("x5", "y5")]) - Sigma

# inspect Psi matrices under mixture trajectory
sim_4_mix$Psi_list[[1]]
sim_4_mix$Psi_list[[2]]
sim_4_mix$Psi_list[[3]]
sim_4_mix$Psi_list[[4]]
sim_4_mix$Psi_list[[5]]

###########################################################################
# VACUUM TEST OF THE STATE VARIANCE EQUATION: MIXTURE TRAJECTORY
###########################################################################

Omega_full_mix <- sim_4_mix$Omega_full

# wave 1
Delta1_mix <- Delta_list_4_mix[[1]]
Psi1_mix <- sim_4_mix$Psi_list[[1]]

Var_W1_mix <- Delta1_mix %*% Omega_full_mix %*% t(Delta1_mix) + Psi1_mix
Var_W1_mix
Var_W1_mix - Sigma

# wave 2
Delta2_mix <- Delta_list_4_mix[[2]]
M1_mix <- sim_4_mix$M_list[[1]]
Psi2_mix <- sim_4_mix$Psi_list[[2]]

Var_W2_mix <-
  Phi %*% Sigma %*% t(Phi) +
  Delta2_mix %*% Omega_full_mix %*% t(Delta2_mix) +
  Phi %*% M1_mix %*% t(Delta2_mix) +
  Delta2_mix %*% t(M1_mix) %*% t(Phi) +
  Psi2_mix

Var_W2_mix
Var_W2_mix - Sigma

# wave 3
Delta3_mix <- Delta_list_4_mix[[3]]
M2_mix <- sim_4_mix$M_list[[2]]
Psi3_mix <- sim_4_mix$Psi_list[[3]]

Var_W3_mix <-
  Phi %*% Sigma %*% t(Phi) +
  Delta3_mix %*% Omega_full_mix %*% t(Delta3_mix) +
  Phi %*% M2_mix %*% t(Delta3_mix) +
  Delta3_mix %*% t(M2_mix) %*% t(Phi) +
  Psi3_mix

Var_W3_mix
Var_W3_mix - Sigma

# wave 4
Delta4_mix <- Delta_list_4_mix[[4]]
M3_mix <- sim_4_mix$M_list[[3]]
Psi4_mix <- sim_4_mix$Psi_list[[4]]

Var_W4_mix <-
  Phi %*% Sigma %*% t(Phi) +
  Delta4_mix %*% Omega_full_mix %*% t(Delta4_mix) +
  Phi %*% M3_mix %*% t(Delta4_mix) +
  Delta4_mix %*% t(M3_mix) %*% t(Phi) +
  Psi4_mix

Var_W4_mix
Var_W4_mix - Sigma

# wave 5
Delta5_mix <- Delta_list_4_mix[[5]]
M4_mix <- sim_4_mix$M_list[[4]]
Psi5_mix <- sim_4_mix$Psi_list[[5]]

Var_W5_mix <-
  Phi %*% Sigma %*% t(Phi) +
  Delta5_mix %*% Omega_full_mix %*% t(Delta5_mix) +
  Phi %*% M4_mix %*% t(Delta5_mix) +
  Delta5_mix %*% t(M4_mix) %*% t(Phi) +
  Psi5_mix

Var_W5_mix
Var_W5_mix - Sigma

###########################################################################
# EMPIRICAL CHECK OF STATIONARY WAVE COVARIANCE FOR SMALLER N
###########################################################################

sim_4_small <- simulate_panel_data(
  N = N_small,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 123
)

cov(sim_4_small$data[, c("x1", "y1")])
cov(sim_4_small$data[, c("x2", "y2")])
cov(sim_4_small$data[, c("x3", "y3")])
cov(sim_4_small$data[, c("x4", "y4")])
cov(sim_4_small$data[, c("x5", "y5")])

# with smaller N the match will be noisier, but still reasonably close
cov(sim_4_small$data[, c("x1", "y1")]) - Sigma
cov(sim_4_small$data[, c("x2", "y2")]) - Sigma
cov(sim_4_small$data[, c("x3", "y3")]) - Sigma
cov(sim_4_small$data[, c("x4", "y4")]) - Sigma
cov(sim_4_small$data[, c("x5", "y5")]) - Sigma

###########################################################################
# REPRODUCIBILITY CHECK
###########################################################################

sim_4_a <- simulate_panel_data(
  N = 5000,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 999
)

sim_4_b <- simulate_panel_data(
  N = 5000,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 999
)

sim_4_c <- simulate_panel_data(
  N = 5000,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 1000
)

identical(sim_4_a$data, sim_4_b$data)
identical(sim_4_a$data, sim_4_c$data)

###########################################################################
# ERROR CHECKS
###########################################################################

# N invalid
try(simulate_panel_data(
  N = 0,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma
))

# T invalid
try(simulate_panel_data(
  N = 100,
  T = 0,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma
))

# Phi wrong dimension
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = diag(3),
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = Sigma
))

# Sigma wrong dimension
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = diag(3)
))

# Sigma not symmetric
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = matrix(c(1, 0.2, 0.1, 1), 2, 2, byrow = TRUE)
))

# Sigma diagonal not equal to 1
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = Omega11_4,
  Sigma = matrix(c(2, 0.3, 0.3, 1), 2, 2, byrow = TRUE)
))

# Omega11 not symmetric
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = matrix(c(1, 0.2, 0.1, 1), 2, 2, byrow = TRUE),
  Sigma = Sigma
))

# Omega11 diagonal not all 1
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const,
  Omega11 = matrix(c(2, 0.3, 0.3, 1), 2, 2, byrow = TRUE),
  Sigma = Sigma
))

# Delta_list wrong length
try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_4_const[1:3],
  Omega11 = Omega11_4,
  Sigma = Sigma
))

# Delta matrices without column names
Delta_list_no_names <- Delta_list_4_const
colnames(Delta_list_no_names[[1]]) <- NULL

try(simulate_panel_data(
  N = 100,
  T = T_panel,
  Phi = Phi,
  Delta_list = Delta_list_no_names,
  Omega11 = Omega11_4,
  Sigma = Sigma
))

###########################################################################
# STRESS CHECK: FEATURE-COVARIANCE LOGIC IN A VACUUM
###########################################################################

# k = 2, independent
sim_cov_2 <- simulate_panel_data(
  N = N_big,
  T = 1,
  Phi = Phi,
  Delta_list = list(Delta1_2),
  Omega11 = Omega11_2,
  Sigma = Sigma,
  seed = 123
)

colnames(sim_cov_2$confounders)
cov(sim_cov_2$confounders)
sim_cov_2$Omega_full
cov(sim_cov_2$confounders) - sim_cov_2$Omega_full

# k = 3, independent
sim_cov_3 <- simulate_panel_data(
  N = N_big,
  T = 1,
  Phi = Phi,
  Delta_list = list(Delta1_3),
  Omega11 = Omega11_3,
  Sigma = Sigma,
  seed = 123
)

colnames(sim_cov_3$confounders)
cov(sim_cov_3$confounders)
sim_cov_3$Omega_full
cov(sim_cov_3$confounders) - sim_cov_3$Omega_full

# k = 4, correlated
sim_cov_4 <- simulate_panel_data(
  N = N_big,
  T = 1,
  Phi = Phi,
  Delta_list = list(Delta1_4),
  Omega11 = Omega11_4,
  Sigma = Sigma,
  seed = 123
)

colnames(sim_cov_4$confounders)
cov(sim_cov_4$confounders)
sim_cov_4$Omega_full
cov(sim_cov_4$confounders) - sim_cov_4$Omega_full