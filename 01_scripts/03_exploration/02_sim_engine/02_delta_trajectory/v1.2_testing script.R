# goal is to stress test all parts of the delta trajectory generation process
# to make sure that everything is working as expected.

###########################################################################
# DELTA TRAJECTORIES
###########################################################################

library(here)
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "02_delta_trajectory", "v1_delta_trajectory.R"))

###########################################################################
# SETUP: BUILD BASELINE OBJECTS
###########################################################################

# independent confounders, k = 1
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
Omega_full_1 <- base_1$Omega

Delta1_1
Omega_full_1


# independent confounders, k = 2
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
Omega_full_2 <- base_2$Omega

Delta1_2
Omega_full_2


# independent confounders, k = 3
Omega11_1 <- diag(3)

base_3 <- sample_delta_1(
  k = 3,
  Omega11 = Omega11_1,
  R2_total = 0.15,
  R2_interaction = 0.30,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

Delta1_3 <- base_3$Delta
Omega_full_3 <- base_3$Omega

Delta1_3
Omega_full_3


# correlated confounders, k = 4
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
Omega_full_4 <- base_4$Omega

Delta1_4
Omega_full_4

###########################################################################
# SPECIAL CASE: k = 1
###########################################################################

Delta1_1
colnames(Delta1_1)
Omega_full_1

# there should only be one coefficient
ncol(Delta1_1)
colnames(Delta1_1)

# baseline implied variances
as.numeric(t(Delta1_1["X", ]) %*% Omega_full_1 %*% Delta1_1["X", ])
as.numeric(t(Delta1_1["Y", ]) %*% Omega_full_1 %*% Delta1_1["Y", ])

# constant trajectory
traj_const_1k <- generate_Delta_constant(
  Delta1 = Delta1_1,
  T = 4
)

traj_const_1k

identical(traj_const_1k[[1]], Delta1_1)
identical(traj_const_1k[[2]], Delta1_1)
identical(traj_const_1k[[3]], Delta1_1)
identical(traj_const_1k[[4]], Delta1_1)

as.numeric(t(traj_const_1k[[1]]["X", ]) %*% Omega_full_1 %*% traj_const_1k[[1]]["X", ])
as.numeric(t(traj_const_1k[[4]]["Y", ]) %*% Omega_full_1 %*% traj_const_1k[[4]]["Y", ])

# stepwise trajectory
traj_step_1k <- generate_Delta_stepwise(
  Delta1 = Delta1_1,
  T = 5,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.40
)

traj_step_1k

identical(traj_step_1k[[1]], Delta1_1)
identical(traj_step_1k[[2]], Delta1_1)

scale_factor_1k <- sqrt(0.40 / 0.15)
scale_factor_1k

all.equal(traj_step_1k[[3]], Delta1_1 * scale_factor_1k)
all.equal(traj_step_1k[[4]], Delta1_1 * scale_factor_1k)
all.equal(traj_step_1k[[5]], Delta1_1 * scale_factor_1k)

abs(as.numeric(t(traj_step_1k[[2]]["X", ]) %*% Omega_full_1 %*% traj_step_1k[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_step_1k[[2]]["Y", ]) %*% Omega_full_1 %*% traj_step_1k[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_step_1k[[4]]["X", ]) %*% Omega_full_1 %*% traj_step_1k[[4]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_step_1k[[4]]["Y", ]) %*% Omega_full_1 %*% traj_step_1k[[4]]["Y", ]) - 0.40) < 1e-8

# mixture trajectory
# here the function should still work, but only through the linear block
traj_mix_1k <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_1,
  T = 5,
  Omega = Omega_full_1,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.50,
  lambda_NL = 0.50,
  seed = 123
)

traj_mix_1k

identical(traj_mix_1k[[1]], Delta1_1)
identical(traj_mix_1k[[2]], Delta1_1)

all.equal(traj_mix_1k[[3]], Delta1_1)
all.equal(traj_mix_1k[[4]], Delta1_1)
all.equal(traj_mix_1k[[5]], Delta1_1)

abs(as.numeric(t(traj_mix_1k[[2]]["X", ]) %*% Omega_full_1 %*% traj_mix_1k[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_mix_1k[[2]]["Y", ]) %*% Omega_full_1 %*% traj_mix_1k[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_mix_1k[[4]]["X", ]) %*% Omega_full_1 %*% traj_mix_1k[[4]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_mix_1k[[4]]["Y", ]) %*% Omega_full_1 %*% traj_mix_1k[[4]]["Y", ]) - 0.40) < 1e-8

###########################################################################
# SPECIAL CASE: k = 2
###########################################################################

Delta1_2
colnames(Delta1_2)
Omega_full_2

# there should be two linear terms and one 2-way interaction
ncol(Delta1_2)
colnames(Delta1_2)

# baseline implied variances
as.numeric(t(Delta1_2["X", ]) %*% Omega_full_2 %*% Delta1_2["X", ])
as.numeric(t(Delta1_2["Y", ]) %*% Omega_full_2 %*% Delta1_2["Y", ])

# constant trajectory
traj_const_2k <- generate_Delta_constant(
  Delta1 = Delta1_2,
  T = 4
)

traj_const_2k

identical(traj_const_2k[[1]], Delta1_2)
identical(traj_const_2k[[2]], Delta1_2)
identical(traj_const_2k[[3]], Delta1_2)
identical(traj_const_2k[[4]], Delta1_2)

as.numeric(t(traj_const_2k[[1]]["X", ]) %*% Omega_full_2 %*% traj_const_2k[[1]]["X", ])
as.numeric(t(traj_const_2k[[4]]["Y", ]) %*% Omega_full_2 %*% traj_const_2k[[4]]["Y", ])

# stepwise trajectory
traj_step_2k <- generate_Delta_stepwise(
  Delta1 = Delta1_2,
  T = 5,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.40
)

traj_step_2k

identical(traj_step_2k[[1]], Delta1_2)
identical(traj_step_2k[[2]], Delta1_2)

scale_factor_2k <- sqrt(0.40 / 0.15)
scale_factor_2k

all.equal(traj_step_2k[[3]], Delta1_2 * scale_factor_2k)
all.equal(traj_step_2k[[4]], Delta1_2 * scale_factor_2k)
all.equal(traj_step_2k[[5]], Delta1_2 * scale_factor_2k)

abs(as.numeric(t(traj_step_2k[[2]]["X", ]) %*% Omega_full_2 %*% traj_step_2k[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_step_2k[[2]]["Y", ]) %*% Omega_full_2 %*% traj_step_2k[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_step_2k[[4]]["X", ]) %*% Omega_full_2 %*% traj_step_2k[[4]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_step_2k[[4]]["Y", ]) %*% Omega_full_2 %*% traj_step_2k[[4]]["Y", ]) - 0.40) < 1e-8

# mixture trajectory
traj_mix_2k <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_2,
  T = 5,
  Omega = Omega_full_2,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.50,
  lambda_NL = 0.50,
  seed = 123
)

traj_mix_2k

identical(traj_mix_2k[[1]], Delta1_2)
identical(traj_mix_2k[[2]], Delta1_2)

all.equal(traj_mix_2k[[3]], Delta1_2)
all.equal(traj_mix_2k[[4]], Delta1_2)
all.equal(traj_mix_2k[[5]], Delta1_2)

abs(as.numeric(t(traj_mix_2k[[2]]["X", ]) %*% Omega_full_2 %*% traj_mix_2k[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_mix_2k[[2]]["Y", ]) %*% Omega_full_2 %*% traj_mix_2k[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_mix_2k[[4]]["X", ]) %*% Omega_full_2 %*% traj_mix_2k[[4]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_mix_2k[[4]]["Y", ]) %*% Omega_full_2 %*% traj_mix_2k[[4]]["Y", ]) - 0.40) < 1e-8

###########################################################################
# EXPLICIT LINEAR / NONLINEAR CHECK FOR k = 2 MIXTURE CASE
###########################################################################

colnames(Delta1_2)

n_colons_2 <- nchar(gsub("[^:]", "", colnames(Delta1_2)))
n_colons_2

idx_L_2 <- which(n_colons_2 == 0)
idx_NL_2 <- which(n_colons_2 >= 1)

idx_L_2
idx_NL_2

d_old_X_2 <- traj_mix_2k[[1]]["X", ]
d_new_X_2 <- traj_mix_2k[[3]]["X", ]

d_old_X_2
d_new_X_2

d_old_X_2_L  <- d_old_X_2[idx_L_2]
d_old_X_2_NL <- d_old_X_2[idx_NL_2]

d_new_X_2_L  <- d_new_X_2[idx_L_2]
d_new_X_2_NL <- d_new_X_2[idx_NL_2]

Omega_L_2   <- Omega_full_2[idx_L_2, idx_L_2, drop = FALSE]
Omega_NL_2  <- Omega_full_2[idx_NL_2, idx_NL_2, drop = FALSE]
Omega_LNL_2 <- Omega_full_2[idx_L_2, idx_NL_2, drop = FALSE]

V_old_X_2_L <- as.numeric(t(d_old_X_2_L) %*% Omega_L_2 %*% d_old_X_2_L)
V_old_X_2_L

V_old_X_2_NL <- as.numeric(t(d_old_X_2_NL) %*% Omega_NL_2 %*% d_old_X_2_NL)
V_old_X_2_NL

V_old_X_2_cross_half <- as.numeric(t(d_old_X_2_L) %*% Omega_LNL_2 %*% d_old_X_2_NL)
V_old_X_2_cross_half

V_old_X_2_L_eff <- V_old_X_2_L + V_old_X_2_cross_half
V_old_X_2_NL_eff <- V_old_X_2_NL + V_old_X_2_cross_half
V_old_X_2_total <- V_old_X_2_L_eff + V_old_X_2_NL_eff

V_old_X_2_L_eff
V_old_X_2_NL_eff
V_old_X_2_total

V_new_X_2_L <- as.numeric(t(d_new_X_2_L) %*% Omega_L_2 %*% d_new_X_2_L)
V_new_X_2_L

V_new_X_2_NL <- as.numeric(t(d_new_X_2_NL) %*% Omega_NL_2 %*% d_new_X_2_NL)
V_new_X_2_NL

V_new_X_2_cross_half <- as.numeric(t(d_new_X_2_L) %*% Omega_LNL_2 %*% d_new_X_2_NL)
V_new_X_2_cross_half

V_new_X_2_L_eff <- V_new_X_2_L + V_new_X_2_cross_half
V_new_X_2_NL_eff <- V_new_X_2_NL + V_new_X_2_cross_half
V_new_X_2_total <- V_new_X_2_L_eff + V_new_X_2_NL_eff

V_new_X_2_L_eff
V_new_X_2_NL_eff
V_new_X_2_total

###########################################################################
# FUNCTION (1): CONSTANT DELTA TRAJECTORY
###########################################################################

traj_const <- generate_Delta_constant(
  Delta1 = Delta1_3,
  T = 5
)

traj_const

# check names
names(traj_const)

# check that each entry is exactly Delta1
identical(traj_const[[1]], Delta1_3)
identical(traj_const[[2]], Delta1_3)
identical(traj_const[[3]], Delta1_3)
identical(traj_const[[4]], Delta1_3)
identical(traj_const[[5]], Delta1_3)

# check row-wise implied variances at first time point
as.numeric(t(traj_const[[1]]["X", ]) %*% Omega_full_3 %*% traj_const[[1]]["X", ])
as.numeric(t(traj_const[[1]]["Y", ]) %*% Omega_full_3 %*% traj_const[[1]]["Y", ])

# check row-wise implied variances at last time point
as.numeric(t(traj_const[[5]]["X", ]) %*% Omega_full_3 %*% traj_const[[5]]["X", ])
as.numeric(t(traj_const[[5]]["Y", ]) %*% Omega_full_3 %*% traj_const[[5]]["Y", ])

###########################################################################
# FUNCTION (2): STEPWISE DELTA TRAJECTORY
###########################################################################

traj_step <- generate_Delta_stepwise(
  Delta1 = Delta1_3,
  T = 6,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

traj_step

# inspect before and after the step
traj_step[[1]]
traj_step[[3]]
traj_step[[4]]
traj_step[[6]]

# before step: should be exactly baseline
identical(traj_step[[1]], Delta1_3)
identical(traj_step[[2]], Delta1_3)
identical(traj_step[[3]], Delta1_3)

# compute expected scaling factor
scale_factor <- sqrt(0.40 / 0.15)
scale_factor

# after step: should be exactly baseline times scaling factor
all.equal(traj_step[[4]], Delta1_3 * scale_factor)
all.equal(traj_step[[5]], Delta1_3 * scale_factor)
all.equal(traj_step[[6]], Delta1_3 * scale_factor)

# check implied row variances before step
as.numeric(t(traj_step[[2]]["X", ]) %*% Omega_full_3 %*% traj_step[[2]]["X", ])
as.numeric(t(traj_step[[2]]["Y", ]) %*% Omega_full_3 %*% traj_step[[2]]["Y", ])

# check implied row variances after step
as.numeric(t(traj_step[[5]]["X", ]) %*% Omega_full_3 %*% traj_step[[5]]["X", ])
as.numeric(t(traj_step[[5]]["Y", ]) %*% Omega_full_3 %*% traj_step[[5]]["Y", ])

# check directly against targets
abs(as.numeric(t(traj_step[[2]]["X", ]) %*% Omega_full_3 %*% traj_step[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_step[[2]]["Y", ]) %*% Omega_full_3 %*% traj_step[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_step[[5]]["X", ]) %*% Omega_full_3 %*% traj_step[[5]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_step[[5]]["Y", ]) %*% Omega_full_3 %*% traj_step[[5]]["Y", ]) - 0.40) < 1e-8

###########################################################################
# FUNCTION (3): STEPWISE MIXTURE DELTA TRAJECTORY
###########################################################################

traj_mix <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_4,
  T = 6,
  Omega = Omega_full_4,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.50,
  lambda_NL = 0.50,
  seed = 123
)

traj_mix

# inspect matrices before and after the step
traj_mix[[1]]
traj_mix[[3]]
traj_mix[[4]]
traj_mix[[6]]

# before step: should be exactly baseline
identical(traj_mix[[1]], Delta1_4)
identical(traj_mix[[2]], Delta1_4)
identical(traj_mix[[3]], Delta1_4)

# after step: generally should not be identical to baseline
all.equal(traj_mix[[4]], Delta1_4)
all.equal(traj_mix[[5]], Delta1_4)
all.equal(traj_mix[[6]], Delta1_4)

# check implied row variances before step
as.numeric(t(traj_mix[[2]]["X", ]) %*% Omega_full_4 %*% traj_mix[[2]]["X", ])
as.numeric(t(traj_mix[[2]]["Y", ]) %*% Omega_full_4 %*% traj_mix[[2]]["Y", ])

# check implied row variances after step
as.numeric(t(traj_mix[[5]]["X", ]) %*% Omega_full_4 %*% traj_mix[[5]]["X", ])
as.numeric(t(traj_mix[[5]]["Y", ]) %*% Omega_full_4 %*% traj_mix[[5]]["Y", ])

# check directly against targets
abs(as.numeric(t(traj_mix[[2]]["X", ]) %*% Omega_full_4 %*% traj_mix[[2]]["X", ]) - 0.15) < 1e-8
abs(as.numeric(t(traj_mix[[2]]["Y", ]) %*% Omega_full_4 %*% traj_mix[[2]]["Y", ]) - 0.15) < 1e-8

abs(as.numeric(t(traj_mix[[5]]["X", ]) %*% Omega_full_4 %*% traj_mix[[5]]["X", ]) - 0.40) < 1e-8
abs(as.numeric(t(traj_mix[[5]]["Y", ]) %*% Omega_full_4 %*% traj_mix[[5]]["Y", ]) - 0.40) < 1e-8

###########################################################################
# EXPLICIT CHECK OF LINEAR VS NONLINEAR PARTS FOR THE MIXTURE FUNCTION
###########################################################################

# identify which columns are linear and which are nonlinear
colnames(Delta1_4)

n_colons <- nchar(gsub("[^:]", "", colnames(Delta1_4)))
n_colons

idx_L <- which(n_colons == 0)
idx_NL <- which(n_colons >= 1)

idx_L
idx_NL

# inspect X row before the step
d_old_X <- traj_mix[[1]]["X", ]
d_new_X <- traj_mix[[4]]["X", ]

d_old_X
d_new_X

# split into linear and nonlinear parts
d_old_X_L  <- d_old_X[idx_L]
d_old_X_NL <- d_old_X[idx_NL]

d_new_X_L  <- d_new_X[idx_L]
d_new_X_NL <- d_new_X[idx_NL]

# build covariance blocks
Omega_L   <- Omega_full_4[idx_L, idx_L, drop = FALSE]
Omega_NL  <- Omega_full_4[idx_NL, idx_NL, drop = FALSE]
Omega_LNL <- Omega_full_4[idx_L, idx_NL, drop = FALSE]

# OLD: linear part
V_old_X_L <- as.numeric(t(d_old_X_L) %*% Omega_L %*% d_old_X_L)
V_old_X_L

# OLD: nonlinear part
V_old_X_NL <- as.numeric(t(d_old_X_NL) %*% Omega_NL %*% d_old_X_NL)
V_old_X_NL

# OLD: half of cross-part
V_old_X_cross_half <- as.numeric(t(d_old_X_L) %*% Omega_LNL %*% d_old_X_NL)
V_old_X_cross_half

# OLD effective decomposition
V_old_X_L_eff  <- V_old_X_L  + V_old_X_cross_half
V_old_X_NL_eff <- V_old_X_NL + V_old_X_cross_half
V_old_X_total  <- V_old_X_L_eff + V_old_X_NL_eff

V_old_X_L_eff
V_old_X_NL_eff
V_old_X_total

# NEW: linear part
V_new_X_L <- as.numeric(t(d_new_X_L) %*% Omega_L %*% d_new_X_L)
V_new_X_L

# NEW: nonlinear part
V_new_X_NL <- as.numeric(t(d_new_X_NL) %*% Omega_NL %*% d_new_X_NL)
V_new_X_NL

# NEW: half of cross-part
V_new_X_cross_half <- as.numeric(t(d_new_X_L) %*% Omega_LNL %*% d_new_X_NL)
V_new_X_cross_half

# NEW effective decomposition
V_new_X_L_eff  <- V_new_X_L  + V_new_X_cross_half
V_new_X_NL_eff <- V_new_X_NL + V_new_X_cross_half
V_new_X_total  <- V_new_X_L_eff + V_new_X_NL_eff

V_new_X_L_eff
V_new_X_NL_eff
V_new_X_total

###########################################################################
# EDGE CASES
###########################################################################

# T = 1 for constant function
traj_const_1 <- generate_Delta_constant(
  Delta1 = Delta1_3,
  T = 1
)

traj_const_1
identical(traj_const_1[[1]], Delta1_3)

# step_at = 1 for stepwise function
traj_step_early <- generate_Delta_stepwise(
  Delta1 = Delta1_3,
  T = 5,
  step_at = 1,
  old_R2 = 0.15,
  new_R2 = 0.30
)

traj_step_early[[1]]
traj_step_early[[5]]

abs(as.numeric(t(traj_step_early[[1]]["X", ]) %*% Omega_full_3 %*% traj_step_early[[1]]["X", ]) - 0.30) < 1e-8
abs(as.numeric(t(traj_step_early[[1]]["Y", ]) %*% Omega_full_3 %*% traj_step_early[[1]]["Y", ]) - 0.30) < 1e-8

# step_at = T for stepwise function
traj_step_late <- generate_Delta_stepwise(
  Delta1 = Delta1_3,
  T = 5,
  step_at = 5,
  old_R2 = 0.15,
  new_R2 = 0.30
)

traj_step_late[[1]]
traj_step_late[[4]]
traj_step_late[[5]]

identical(traj_step_late[[1]], Delta1_3)
identical(traj_step_late[[4]], Delta1_3)

abs(as.numeric(t(traj_step_late[[5]]["X", ]) %*% Omega_full_3 %*% traj_step_late[[5]]["X", ]) - 0.30) < 1e-8
abs(as.numeric(t(traj_step_late[[5]]["Y", ]) %*% Omega_full_3 %*% traj_step_late[[5]]["Y", ]) - 0.30) < 1e-8

# mixture with lambda = 0
traj_mix_zero <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_4,
  T = 5,
  Omega = Omega_full_4,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.30,
  lambda_L = 0,
  lambda_NL = 0,
  seed = 123
)

traj_mix_zero[[2]]
traj_mix_zero[[3]]

abs(as.numeric(t(traj_mix_zero[[3]]["X", ]) %*% Omega_full_4 %*% traj_mix_zero[[3]]["X", ]) - 0.30) < 1e-8
abs(as.numeric(t(traj_mix_zero[[3]]["Y", ]) %*% Omega_full_4 %*% traj_mix_zero[[3]]["Y", ]) - 0.30) < 1e-8

# mixture with lambda = 1
traj_mix_one <- generate_Delta_stepwise_mixture(
  Delta1 = Delta1_4,
  T = 5,
  Omega = Omega_full_4,
  step_at = 3,
  old_R2 = 0.15,
  new_R2 = 0.30,
  lambda_L = 1,
  lambda_NL = 1,
  seed = 123
)

traj_mix_one[[2]]
traj_mix_one[[3]]

abs(as.numeric(t(traj_mix_one[[3]]["X", ]) %*% Omega_full_4 %*% traj_mix_one[[3]]["X", ]) - 0.30) < 1e-8
abs(as.numeric(t(traj_mix_one[[3]]["Y", ]) %*% Omega_full_4 %*% traj_mix_one[[3]]["Y", ]) - 0.30) < 1e-8
