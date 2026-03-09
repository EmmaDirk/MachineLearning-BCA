# this script serves to test the trajectory functions.
#
#
# Recall that the trajectory generators:
# - take a baseline sampled Delta-matrix D1,
# - generate a list of D-matrices across T time points,
# - either keep the confounder effects constant,
# - or change them stepwise by a common rescaling,
# - or change them stepwise by a convex mixture plus exact rescaling,
# - preserve matrix structure and naming across time,
# - and error when logical input constraints are violated.
#
# We therefore test:
#
# 1) That the trajectory generators error under invalid inputs:
#    a) D1 not a matrix.
#    b) T invalid.
#    c) step_at invalid for the stepwise generators.
#    d) old_R2 / new_R2 invalid for the stepwise generators.
#    e) Omega missing or invalid for the mixture generator.
#    f) lambda_L / lambda_NL invalid for the mixture generator.
#    g) D1 missing column names for the mixture generator.
#
# 2) That the returned objects have the correct structure:
#    a) A list.
#    b) Of length T.
#    c) Named t1, ..., tT.
#    d) Each element is a matrix with the same dimensions as D1.
#    e) Row names and column names of D1 are preserved.
#
# 3) That the scenarios behave as intended:
#    a) Constant scenario returns D1 at every time point.
#    b) Stepwise scenario returns D1 before the step.
#    c) Stepwise scenario returns D1 * sqrt(new_R2 / old_R2) after the step.
#    d) Stepwise mixture returns D1 before the step.
#    e) Stepwise mixture returns one fixed post-step matrix after the step.
#    f) Stepwise mixture generally differs from simple common rescaling.
#
# 4) That the mathematical targets are satisfied:
#    a) Stepwise scenario hits the intended post-step R2.
#    b) Stepwise mixture hits old_R2 at baseline.
#    c) Stepwise mixture hits new_R2 after the step.
#    d) When nonlinear terms are present, the linear / nonlinear variance shares are preserved.
#
# 5) That the generators are reproducible when appropriate:
#    a) Constant generator is deterministic.
#    b) Stepwise generator is deterministic.
#    c) Mixture generator is reproducible when the seed is fixed.
#
# ------------------------------------------------------------------------------------------------------------------

# libraries
library(here)
library(ggplot2)
library(viridis)

# source scripts
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "02_sim_engine", "02_delta_trajectory", "v1_delta_trajectory.R"))

# reproducibility
set.seed(123)

# ------------------------------------------------------------------------------------------------------------------
# helpers used in the script
# ------------------------------------------------------------------------------------------------------------------

# helper to compute quadratic-form variance
qvar <- function(d, Om) {
  as.numeric(t(d) %*% Om %*% d)
}

# helper to count colons in coefficient names
count_colons <- function(x) {
  nchar(gsub("[^:]", "", x))
}

# helper to convert a D_list into long format for plotting
D_list_to_long <- function(D_list, scenario_name) {

  # build a list of data frames
  out <- vector("list", length(D_list))

  # fill the list
  for (t in seq_along(D_list)) {

    # extract D_t
    D_t <- D_list[[t]]

    # convert to a data frame
    df_t <- data.frame(

      # where we have the collumns:
      time = t,
      outcome = rep(rownames(D_t), each = ncol(D_t)),
      coefficient = rep(colnames(D_t), times = nrow(D_t)),
      delta = as.vector(t(D_t)),
      scenario = scenario_name,
      stringsAsFactors = FALSE
    )

    # add to the list
    out[[t]] <- df_t
  }

  # bind the data frames
  do.call(rbind, out)
}

# ------------------------------------------------------------------------------------------------------------------
# create an example D1 and Omega used throughout the checks
# ------------------------------------------------------------------------------------------------------------------

# create simple D1
out1 <- sample_delta_1(
  k = 3,
  Sigma = diag(3),
  R2_total = 0.15,
  R2_interaction = 0.30,
  include_2way = TRUE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

# extract elements
D1    <- out1$Delta
Omega <- out1$Omega

# number of time points
T <- 5

# identify linear and nonlinear coefficients
n_colons <- count_colons(colnames(D1))
idx_L  <- which(n_colons == 0)
idx_NL <- which(n_colons >= 1)

# ------------------------------------------------------------------------------------------------------------------
# 0: visualize the delta trajectories under the 3 scenarios
# ------------------------------------------------------------------------------------------------------------------

# constant
D_const <- generate_D_constant(
  D1 = D1,
  T  = T
)

# stepwise
D_step <- generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# stepwise mixture
D_mix <- generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.3,
  lambda_NL = 0.3,

  # important: this seed must differ from the seed used above for sampling D1
  seed = 1234
)

# look at the results
# visually inspect
D_const
D_step
D_mix

# build one combined plotting data set
plot_dat <- rbind(
  D_list_to_long(D_const, "constant"),
  D_list_to_long(D_step,  "stepwise"),
  D_list_to_long(D_mix,   "stepwise mixture")
)

# inspect
plot_dat

# make scenario order explicit
plot_dat$scenario <- factor(
  plot_dat$scenario,
  levels = c("constant", "stepwise", "stepwise mixture")
)

# ensure coefficient colors are stable
plot_dat$coefficient <- factor(
  plot_dat$coefficient,
  levels = unique(plot_dat$coefficient)
)

# plot
ggplot(plot_dat, aes(x = time, y = delta, group = coefficient, color = coefficient)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(outcome ~ scenario, scales = "free_y") +
  labs(
    x = "Time",
    y = expression(delta),
    color = "Coefficient",
    title = "Delta trajectories under three scenarios"
  ) +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------------------------------------------
# 1: test that the trajectory generators error under invalid inputs
# ------------------------------------------------------------------------------------------------------------------

# a1) generate_D_constant(): D1 not a matrix
# SHOULD ERROR
generate_D_constant(
  D1 = 1:4,
  T = 5
)

# a2) generate_D_stepwise(): D1 not a matrix
# SHOULD ERROR
generate_D_stepwise(
  D1 = 1:4,
  T = 5,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# a3) generate_D_stepwise_mixture(): D1 not a matrix
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = 1:4,
  T = 5,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# b1) generate_D_constant(): T invalid
# SHOULD ERROR
generate_D_constant(
  D1 = D1,
  T = 0
)

# b2) generate_D_stepwise(): T invalid
# SHOULD ERROR
generate_D_stepwise(
  D1 = D1,
  T = 0,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# b3) generate_D_stepwise_mixture(): T invalid
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = 0,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# c1) generate_D_stepwise(): step_at outside 1..T
# SHOULD ERROR
generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = T + 1,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# c2) generate_D_stepwise_mixture(): step_at outside 1..T
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = T + 1,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# d1) generate_D_stepwise(): old_R2 invalid
# SHOULD ERROR
generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0,
  new_R2 = 0.40
)

# d2) generate_D_stepwise(): new_R2 invalid
# SHOULD ERROR
generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 2
)

# d3) generate_D_stepwise_mixture(): old_R2 invalid
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0,
  new_R2 = 0.40
)

# d4) generate_D_stepwise_mixture(): new_R2 invalid
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 2
)

# e1) generate_D_stepwise_mixture(): Omega missing
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# e2) generate_D_stepwise_mixture(): Omega not a matrix
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = 1:4,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# e3) generate_D_stepwise_mixture(): Omega not symmetric
# SHOULD ERROR
Omega_nonsym <- Omega
Omega_nonsym[1, 2] <- Omega_nonsym[1, 2] + 0.1

generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega_nonsym,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# e4) generate_D_stepwise_mixture(): D1 and Omega not conformable
# SHOULD ERROR
Omega_small <- Omega[-1, -1]

generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega_small,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# f1) generate_D_stepwise_mixture(): lambda_L invalid
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 2,
  lambda_NL = 0.3
)

# f2) generate_D_stepwise_mixture(): lambda_NL invalid
# SHOULD ERROR
generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.3,
  lambda_NL = -1
)

# g) generate_D_stepwise_mixture(): D1 has no column names
# SHOULD ERROR
D1_noname <- D1
colnames(D1_noname) <- NULL

generate_D_stepwise_mixture(
  D1 = D1_noname,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

# ------------------------------------------------------------------------------------------------------------------
# 2: test that the returned objects have the correct structure
# ------------------------------------------------------------------------------------------------------------------

# visually inspect the 3 outputs
# should be:
# a) a list.
# b) of length T.
# c) named t1, ..., tT.
# d) each element a matrix with same dimensions as D1.
# e) row names and column names preserved.
D_const
D_step
D_mix

# inspect structure
str(D_const)
str(D_step)
str(D_mix)

# basic checks for constant scenario
# SHOULD BE TRUE
is.list(D_const)
length(D_const) == T
identical(names(D_const), paste0("t", 1:T))
is.matrix(D_const[[1]])
identical(dim(D_const[[1]]), dim(D1))
identical(rownames(D_const[[1]]), rownames(D1))
identical(colnames(D_const[[1]]), colnames(D1))

# basic checks for stepwise scenario
# SHOULD BE TRUE
is.list(D_step)
length(D_step) == T
identical(names(D_step), paste0("t", 1:T))
is.matrix(D_step[[1]])
identical(dim(D_step[[1]]), dim(D1))
identical(rownames(D_step[[1]]), rownames(D1))
identical(colnames(D_step[[1]]), colnames(D1))

# basic checks for stepwise mixture scenario
# SHOULD BE TRUE
is.list(D_mix)
length(D_mix) == T
identical(names(D_mix), paste0("t", 1:T))
is.matrix(D_mix[[1]])
identical(dim(D_mix[[1]]), dim(D1))
identical(rownames(D_mix[[1]]), rownames(D1))
identical(colnames(D_mix[[1]]), colnames(D1))

# ------------------------------------------------------------------------------------------------------------------
# 3: test that the scenarios behave as intended
# ------------------------------------------------------------------------------------------------------------------

# a) constant scenario returns D1 at every time point
# SHOULD BE TRUE
all(vapply(D_const, function(x) isTRUE(all.equal(x, D1)), logical(1)))

# b) stepwise scenario returns D1 before the step
# SHOULD BE TRUE
isTRUE(all.equal(D_step[[1]], D1))
isTRUE(all.equal(D_step[[2]], D1))
isTRUE(all.equal(D_step[[3]], D1))

# c) stepwise scenario returns D1 * sqrt(new_R2 / old_R2) after the step
# SHOULD BE TRUE
scale_factor <- sqrt(0.40 / 0.15)
D_target <- D1 * scale_factor

isTRUE(all.equal(D_step[[4]], D_target))
isTRUE(all.equal(D_step[[5]], D_target))

# d) stepwise mixture returns D1 before the step
# SHOULD BE TRUE
isTRUE(all.equal(D_mix[[1]], D1))
isTRUE(all.equal(D_mix[[2]], D1))
isTRUE(all.equal(D_mix[[3]], D1))

# e) stepwise mixture returns one fixed post-step matrix after the step
# SHOULD BE TRUE
isTRUE(all.equal(D_mix[[4]], D_mix[[5]]))

# f) stepwise mixture generally differs from simple common rescaling
# SHOULD BE FALSE for all.equal(), hence !all.equal(...) should be TRUE
!isTRUE(all.equal(D_mix[[4]], D_target))

# g) stepwise scenario leaves everything unchanged when old_R2 = new_R2
# SHOULD BE TRUE
D_step_same <- generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.15
)

all(vapply(D_step_same, function(x) isTRUE(all.equal(x, D1)), logical(1)))

# h) stepwise scenario preserves coefficient ratios exactly
# SHOULD BE TRUE
nz <- abs(D1) > 1e-12
all(abs(D_step[[4]][nz] / D1[nz] - scale_factor) < 1e-10)

# ------------------------------------------------------------------------------------------------------------------
# 4: test that the mathematical targets are satisfied
# ------------------------------------------------------------------------------------------------------------------

# a) stepwise scenario hits the intended post-step R2
# since all coefficients are scaled by sqrt(new_R2 / old_R2),
# the quadratic-form variance should scale from old_R2 to new_R2
# SHOULD BE TRUE
abs(qvar(D_step[[4]]["X", ], Omega) - 0.40) < 1e-8
abs(qvar(D_step[[4]]["Y", ], Omega) - 0.40) < 1e-8

# b) stepwise mixture hits old_R2 at baseline
# SHOULD BE TRUE
abs(qvar(D_mix[[1]]["X", ], Omega) - 0.15) < 1e-8
abs(qvar(D_mix[[1]]["Y", ], Omega) - 0.15) < 1e-8

# c) stepwise mixture hits new_R2 after the step
# SHOULD BE TRUE
abs(qvar(D_mix[[4]]["X", ], Omega) - 0.40) < 1e-8
abs(qvar(D_mix[[4]]["Y", ], Omega) - 0.40) < 1e-8

# d) when nonlinear terms are present, the linear / nonlinear variance shares are preserved
# first compute the shares at baseline and after the step
Omega_L  <- Omega[idx_L, idx_L, drop = FALSE]
Omega_NL <- Omega[idx_NL, idx_NL, drop = FALSE]

# baseline shares for X
V_L_old_X  <- qvar(D1["X", idx_L], Omega_L)
V_NL_old_X <- qvar(D1["X", idx_NL], Omega_NL)
share_L_old_X  <- V_L_old_X  / (V_L_old_X + V_NL_old_X)
share_NL_old_X <- V_NL_old_X / (V_L_old_X + V_NL_old_X)

# post-step shares for X
V_L_new_X  <- qvar(D_mix[[4]]["X", idx_L], Omega_L)
V_NL_new_X <- qvar(D_mix[[4]]["X", idx_NL], Omega_NL)
share_L_new_X  <- V_L_new_X  / (V_L_new_X + V_NL_new_X)
share_NL_new_X <- V_NL_new_X / (V_L_new_X + V_NL_new_X)

# baseline shares for Y
V_L_old_Y  <- qvar(D1["Y", idx_L], Omega_L)
V_NL_old_Y <- qvar(D1["Y", idx_NL], Omega_NL)
share_L_old_Y  <- V_L_old_Y  / (V_L_old_Y + V_NL_old_Y)
share_NL_old_Y <- V_NL_old_Y / (V_L_old_Y + V_NL_old_Y)

# post-step shares for Y
V_L_new_Y  <- qvar(D_mix[[4]]["Y", idx_L], Omega_L)
V_NL_new_Y <- qvar(D_mix[[4]]["Y", idx_NL], Omega_NL)
share_L_new_Y  <- V_L_new_Y  / (V_L_new_Y + V_NL_new_Y)
share_NL_new_Y <- V_NL_new_Y / (V_L_new_Y + V_NL_new_Y)

# SHOULD BE TRUE
abs(share_L_old_X  - share_L_new_X)  < 1e-8
abs(share_NL_old_X - share_NL_new_X) < 1e-8
abs(share_L_old_Y  - share_L_new_Y)  < 1e-8
abs(share_NL_old_Y - share_NL_new_Y) < 1e-8

# ------------------------------------------------------------------------------------------------------------------
# 5: test that the generators are reproducible when appropriate
# ------------------------------------------------------------------------------------------------------------------

# a) constant generator is deterministic
# SHOULD BE TRUE
D_const_2 <- generate_D_constant(
  D1 = D1,
  T = T
)

all.equal(D_const, D_const_2)

# b) stepwise generator is deterministic
# SHOULD BE TRUE
D_step_2 <- generate_D_stepwise(
  D1 = D1,
  T = T,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40
)

all.equal(D_step, D_step_2)

# c) mixture generator is reproducible when the seed is fixed
# SHOULD BE TRUE
D_mix_2 <- generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.3,
  lambda_NL = 0.3,
  seed = 1234
)

all.equal(D_mix, D_mix_2)

# d) mixture generator changes when a different seed is used
# SHOULD BE FALSE for all.equal(), hence !all.equal(...) should be TRUE
D_mix_other <- generate_D_stepwise_mixture(
  D1 = D1,
  T = T,
  Omega = Omega,
  step_at = 4,
  old_R2 = 0.15,
  new_R2 = 0.40,
  lambda_L = 0.3,
  lambda_NL = 0.3,
  seed = 999
)

!isTRUE(all.equal(D_mix, D_mix_other))
