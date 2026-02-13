# This script is for trying out MSMs, with their inverse probability weights to deal with confounding. 
# we will simulate the same data as we did for the simulation studies, namely:
# - a system with time-varying but linear effects of baseline confounders
# - a system with time-varying, non-linear effects of baseline confounders
# there are no time-varying confounders in the data-generating mechanism. 
# ----------------------------------------------------------------------------------------------------------------

#################################
# step 0: simulate the data
#################################

# step 0.1: load the required packages
library(here)
library(tidyverse)
library(mvtnorm)

# step 0.2: load the required functions for data simulation
source(here("01_scripts", "03_exploration", "01_models", "01_data_sim", "01_delta_sampler.R"))
source(here("01_scripts", "03_exploration", "01_models", "01_data_sim", "02_delta_trajectory.R"))
source(here("01_scripts", "03_exploration", "01_models", "01_data_sim", "03.1_simulate_panel_data.R"))
source(here("01_scripts", "03_exploration", "01_models", "01_data_sim", "03.2_simulate_panel_data.R"))

# step 0.3: set seed for reproducibility
set.seed(1234)

# step 0.4: sample a detla matrix for both data-sets
# first for the data-set with linear confounders
D1_lin <- sample_delta_1(
  k = 3,                                               # number of confounders
  R2_total = 0.15,                                     # total confounder R2 at time t = 1
  min_abs = 0.001,                                     # minimum absolute value for each delta
  max_abs = 0.50,                                      # maximum absolute value for each delta
  R2_nonlin = 0                                        # fraction of R2_total allocated to non-linear terms
)

# then for the data-set with non-linear confounders
D1_nonlin <- sample_delta_1(
  k = 3,                                               # number of confounders
  R2_total = 0.15,                                     # total confounder R2 at time t = 1
  min_abs = 0.001,                                     # minimum absolute value for each coefficient
  max_abs = 0.50,                                      # maximum absolute value for each coefficient
  R2_nonlin = 0.4                                      # fraction of R2_1 allocated to non-linear terms
)

# step 0.5: make delta trajectories from the sampled delta matrices

# first for the data-set with linear confounder effects
# make a constant delta trajectory over 5 time points
D_list_constant_lin <- generate_D_constant(
  D1 = D1_lin,                                         # initial delta matrix
  T  = 5                                               # number of time points
)

# make a stepwise delta trajectory over 5 time points
D_list_stepwise_lin <- generate_D_stepwise(
  D1      = D1_lin,                                    # initial delta matrix
  T       = 5,                                         # number of time points
  step_at = 4,                                         # time point at which to step
  old_R2  = 0.15,                                      # old R2 before the step
  new_R2  = 0.40                                       # new R2 after the step
)

# then for the data-set with non-linear confounder effects
# make a constant delta trajectory over 5 time points
D_list_constant_nonlin <- generate_D_constant(
  D1 = D1_nonlin,                                      # initial delta matrix
  T = 5                                                # number of time points
)

# make a stepwise delta trajectory over 5 time points
D_list_stepwise_nonlin <- generate_D_stepwise(
  D1      = D1_nonlin,                                 # initial delta matrix
  T       = 5,                                         # number of time points
  step_at = 4,                                         # time point at which to step
  old_R2  = 0.15,                                      # old R2 before the step
  new_R2  = 0.40                                       # new R2 after the step
)

# step 0.6: define the true A matrix (beta = autoregressive, gamma = cross-lagged)
A <- matrix(c(
  0.20, 0.00,                                          # BetaX,  GammaYX
  0.10, 0.20                                           # GammaXY, BetaY
), nrow = 2, byrow = TRUE)

# step 0.7: define the confounder covariance matrix
Psi_lin <- diag(3)                                     # uncorrelated confounders with var=1
Psi_nonlin <- diag(7)                                  # uncorrelated confounders with var=1

# step 0.8: simulate panel data for the 4 scenarios (constant vs stepwise, linear vs non-linear)
# 0.8.1: constant effects, linear confounders
data_constant_lin <- simulate_panel_data(
  D_list    = D_list_constant_lin,                     # list of D matrices
  A         = A,                                       # A matrix
  Psi       = Psi_lin,                                 # Psi matrix
  N         = 1000000,                                 # sample size
  T         = 5,                                       # number of waves
  rho_extra = 0                                      # residual correlation
)

# 0.8.2: stepwise effects, linear confounders
data_stepwise_lin <- simulate_panel_data(
  D_list    = D_list_stepwise_lin,
  A         = A,
  Psi       = Psi_lin,
  N         = 1000000,
  T         = 5,
  rho_extra = 0
)

# 0.8.3: constant effects, non-linear confounders
data_constant_nonlin <- simulate_panel_data(
  D_list    = D_list_constant_nonlin,
  A         = A,
  Psi       = Psi_nonlin,
  N         = 100000,
  T         = 5,
  rho_extra = 0
)

# 0.8.4: stepwise effects, non-linear confounders
data_stepwise_nonlin <- simulate_panel_data(
  D_list    = D_list_stepwise_nonlin,
  A         = A,
  Psi       = Psi_nonlin,
  N         = 100000,
  T         = 5,
  rho_extra = 0
)

#############################################
# step 1: fit a MSM without correction for C
#############################################

# We first need to define what our interest is. One of the advantages of MSMs is that we do not have to specify the
# entire datagenerating system, like we would in a SEM. Instead we need to identify all potential confounders for the 
# outcome of interest, and then adjust for these confounders using IPW. 
#
# We have 5 datapoint. We would like to know how X_4 influences Y_5. We need to adjust for all confounders of this relation.
# In this case, we can block all backdoor paths by conditioning on Y_4, and C, since the effect of previous exposures
# always goes through Y_4 in our data-generating mechanism. 
# In our most parsimonious model, we thus only have to model:
# - Y_4 -> Y_5
# - X_4 -> Y_5
# - C -> Y_5
# and for more robustness, we could model: C -> X_4 too. However, in IPW we often try to model the probability of being treated.
# we cannot use Y_4 here, since Y_4 does not cause X_4. 
# As such we need to model:
# - C -> X_4
# - X_3 -> X_4
# - Y_3 -> X_4
# 
# We here choose to only model the treatment probabilities. 
#
# At first we will act as though we do not observe C. In this case, we will use IPW only with X_3 and Y_3.
# In later stages, we will add C to the model. We will use stabilized weights, but not truncated weights.

# 1.1: time-invariant effects of linear confounders
# ---------------------------------------------------

# 1.1.1: compute weights
data_constant_lin <- as.data.frame(data_constant_lin)

# numerator: f(X_4)
num_cl <- lm(x4 ~ 1, data = data_constant_lin)
mu_num_cl <- predict(num_cl)
sd_num_cl <- sd(resid(num_cl))

# denominator: f(X_4 | X_3, Y_3)
denom_cl <- lm(x4 ~ x3 + y3, data = data_constant_lin)
mu_denom_cl <- predict(denom_cl)
sd_denom_cl <- sd(resid(denom_cl))

# compute stabilized density ratio-weights
w_cl <- dnorm(data_constant_lin$x4, mean = mu_num_cl, sd = sd_num_cl) / dnorm(data_constant_lin$x4, mean = mu_denom_cl, sd = sd_denom_cl)

# 1.1.2: fit the MSM
msm_cl <- lm(y5 ~ x4, data = data_constant_lin, weights = w_cl)
summary(msm_cl)

# 1.2: time-variant effects of linear confounders
# ---------------------------------------------------

# 1.2.1: compute weights
data_stepwise_lin <- as.data.frame(data_stepwise_lin)

# numerator: f(X_4)
num_sl <- lm(x4 ~ 1, data = data_stepwise_lin)
mu_num_sl <- predict(num_sl)
sd_num_sl <- sd(resid(num_sl))

# denominator: f(X_4 | X_3, Y_3)
denom_sl <- lm(x4 ~ x3 + y3, data = data_stepwise_lin)
mu_denom_sl <- predict(denom_sl)
sd_denom_sl <- sd(resid(denom_sl))

# compute stabilized density ratio-weights
w_sl <- dnorm(data_stepwise_lin$x4, mean = mu_num_sl, sd = sd_num_sl) / dnorm(data_stepwise_lin$x4, mean = mu_denom_sl, sd = sd_denom_sl)

# 1.2.2: fit the MSM  
msm_sl <- lm(y5 ~ x4, data = data_stepwise_lin, weights = w_sl)
summary(msm_sl)

# 1.3: time-invariant effects of non-linear confounders
# ---------------------------------------------------

# 1.3.1: compute weights
data_constant_nonlin <- as.data.frame(data_constant_nonlin)

# numerator: f(X_4)
num_cn <- lm(x4 ~ 1, data = data_constant_nonlin)
mu_num_cn <- predict(num_cn)
sd_num_cn <- sd(resid(num_cn))

# denominator: f(X_4 | X_3, Y_3)
denom_cn <- lm(x4 ~ x3 + y3, data = data_constant_nonlin)
mu_denom_cn <- predict(denom_cn)
sd_denom_cn <- sd(resid(denom_cn))

# compute stabilized density ratio-weights
w_cn <- dnorm(data_constant_nonlin$x4, mean = mu_num_cn, sd = sd_num_cn) / dnorm(data_constant_nonlin$x4, mean = mu_denom_cn, sd = sd_denom_cn)

# 1.3.2: fit the MSM
msm_cn <- lm(y5 ~ x4, data = data_constant_nonlin, weights = w_cn)
summary(msm_cn)

# 1.4: time-variant effects of non-linear confounders
# ---------------------------------------------------

# 1.4.1: compute weights
data_stepwise_nonlin <- as.data.frame(data_stepwise_nonlin)

# numerator: f(X_4)
num_sn <- lm(x4 ~ 1, data = data_stepwise_nonlin)
mu_num_sn <- predict(num_sn)
sd_num_sn <- sd(resid(num_sn))

# denominator: f(X_4 | X_3, Y_3)
denom_sn <- lm(x4 ~ x3 + y3, data = data_stepwise_nonlin)
mu_denom_sn <- predict(denom_sn)
sd_denom_sn <- sd(resid(denom_sn))

# compute stabilized density ratio-weights
w_sn <- dnorm(data_stepwise_nonlin$x4, mean = mu_num_sn, sd = sd_num_sn) / dnorm(data_stepwise_nonlin$x4, mean = mu_denom_sn, sd = sd_denom_sn)

# 1.4.2: fit the MSM
msm_sn <- lm(y5 ~ x4, data = data_stepwise_nonlin, weights = w_sn)
summary(msm_sn)

############################################
# step 2: fit a MSM with correction for C
############################################

# 2.1: time-invariant effects of linear confounders
# ---------------------------------------------------

# 2.1.1: compute weights
data_constant_lin <- as.data.frame(data_constant_lin)

# numerator: f(X_4)
num_clC <- lm(x4 ~ 1, data = data_constant_lin)
mu_num_clC <- predict(num_clC)
sd_num_clC <- sd(resid(num_clC))

# denominator: f(X_4 | X_3, Y_3, C)
denom_clC <- lm(x4 ~ x3 + y3 + c1 + c2 + c3, data = data_constant_lin)
mu_denom_clC <- predict(denom_clC)
sd_denom_clC <- sd(resid(denom_clC))

# compute stabilized density ratio-weights
w_clC <- dnorm(data_constant_lin$x4, mean = mu_num_clC, sd = sd_num_clC) / dnorm(data_constant_lin$x4, mean = mu_denom_clC, sd = sd_denom_clC)

# 2.1.2: fit the MSM
msm_clC <- lm(y5 ~ x4, data = data_constant_lin, weights = w_clC)
summary(msm_clC)

# 2.2: time-variant effects of linear confounders
# ---------------------------------------------------

# 2.2.1: compute weights
data_stepwise_lin <- as.data.frame(data_stepwise_lin)

# numerator: f(X_4)
num_slC <- lm(x4 ~ 1, data = data_stepwise_lin)
mu_num_slC <- predict(num_slC)
sd_num_slC <- sd(resid(num_slC))

# denominator: f(X_4 | X_3, Y_3, C)
denom_slC <- lm(x4 ~ x3 + y3 + c1 + c2 + c3, data = data_stepwise_lin)
mu_denom_slC <- predict(denom_slC)
sd_denom_slC <- sd(resid(denom_slC))

# compute stabilized density ratio-weights
w_slC <- dnorm(data_stepwise_lin$x4, mean = mu_num_slC, sd = sd_num_slC) / dnorm(data_stepwise_lin$x4, mean = mu_denom_slC, sd = sd_denom_slC)

# 2.2.2: fit the MSM
msm_slC <- lm(y5 ~ x4, data = data_stepwise_lin, weights = w_slC)
summary(msm_slC)

# 2.3: time-invariant effects of non-linear confounders
# ---------------------------------------------------

# 2.3.1: compute weights
data_constant_nonlin <- as.data.frame(data_constant_nonlin)

# numerator: f(X_4)
num_cnC <- lm(x4 ~ 1, data = data_constant_nonlin)
mu_num_cnC <- predict(num_cnC)
sd_num_cnC <- sd(resid(num_cnC))

# denominator: f(X_4 | X_3, Y_3, C)
denom_cnC <- lm(x4 ~ x3 + y3 + c1 + c2 + c3, data = data_constant_nonlin)
mu_denom_cnC <- predict(denom_cnC)
sd_denom_cnC <- sd(resid(denom_cnC))

# compute stabilized density ratio-weights
w_cnC <- dnorm(data_constant_nonlin$x4, mean = mu_num_cnC, sd = sd_num_cnC) / dnorm(data_constant_nonlin$x4, mean = mu_denom_cnC, sd = sd_denom_cnC)

# 2.3.2: fit the MSM
msm_cnC <- lm(y5 ~ x4, data = data_constant_nonlin, weights = w_cnC)
summary(msm_cnC)

# 2.4: time-variant effects of non-linear confounders
# ---------------------------------------------------

# 2.4.1: compute weights
data_stepwise_nonlin <- as.data.frame(data_stepwise_nonlin)

# numerator: f(X_4)
num_snC <- lm(x4 ~ 1, data = data_stepwise_nonlin)
mu_num_snC <- predict(num_snC)
sd_num_snC <- sd(resid(num_snC))

# denominator: f(X_4 | X_3, Y_3, C)
denom_snC <- lm(x4 ~ x3 + y3 + c1 + c2 + c3, data = data_stepwise_nonlin)
mu_denom_snC <- predict(denom_snC)
sd_denom_snC <- sd(resid(denom_snC))

# compute stabilized density ratio-weights
w_snC <- dnorm(data_stepwise_nonlin$x4, mean = mu_num_snC, sd = sd_num_snC) / dnorm(data_stepwise_nonlin$x4, mean = mu_denom_snC, sd = sd_denom_snC)

# 2.4.2: fit the MSM
msm_snC <- lm(y5 ~ x4, data = data_stepwise_nonlin, weights = w_snC)
summary(msm_snC) 

#####################################################################################
# step 3: modelling the interaction between confounders (only nonlinear confounders)
# here we are applying the correct model
#####################################################################################

# 3.1: time-invariant effects of non-linear confounders
# ---------------------------------------------------

# 3.1.1: compute weights
data_constant_nonlin <- as.data.frame(data_constant_nonlin)

# numerator: f(X_4)
num_cnCI <- lm(x4 ~ 1, data = data_constant_nonlin)
mu_num_cnCI <- predict(num_cnCI)
sd_num_cnCI <- sd(resid(num_cnCI))

# denominator: f(X_4 | X_3, Y_3, C, C^2, C^3)
denom_cnCI <- lm(x4 ~ x3 + y3 + c1 + c2 + c3 + c4 + c5 + c6 + c7, data = data_constant_nonlin)
mu_denom_cnCI <- predict(denom_cnCI)
sd_denom_cnCI <- sd(resid(denom_cnCI))

# compute stabilized density ratio-weights
w_cnCI <- dnorm(data_constant_nonlin$x4, mean = mu_num_cnCI, sd = sd_num_cnCI) / dnorm(data_constant_nonlin$x4, mean = mu_denom_cnCI, sd = sd_denom_cnCI)

# 3.1.2: fit the MSM
msm_cnCI <- lm(y5 ~ x4 + y4, data = data_constant_nonlin, weights = w_cnCI)
summary(msm_cnCI)

# 3.2: time-variant effects of non-linear confounders
# ---------------------------------------------------

# 3.2.1: compute weights
data_stepwise_nonlin <- as.data.frame(data_stepwise_nonlin)

# numerator: f(X_4)
num_snCI <- lm(x4 ~ 1, data = data_stepwise_nonlin)
mu_num_snCI <- predict(num_snCI)
sd_num_snCI <- sd(resid(num_snCI))

# denominator: f(X_4 | X_3, Y_3, C, C^2, C^3)
denom_snCI <- lm(x4 ~ x3 + y3 + c1 + c2 + c3 + c4 + c5 + c6 + c7, data = data_stepwise_nonlin)
mu_denom_snCI <- predict(denom_snCI)
sd_denom_snCI <- sd(resid(denom_snCI))

# compute stabilized density ratio-weights
w_snCI <- dnorm(data_stepwise_nonlin$x4, mean = mu_num_snCI, sd = sd_num_snCI) / dnorm(data_stepwise_nonlin$x4, mean = mu_denom_snCI, sd = sd_denom_snCI)

# 3.2.2: fit the MSM
msm_snCI <- lm(y5 ~ x4 + y4, data = data_stepwise_nonlin, weights = w_snCI)
summary(msm_snCI)

# Quick overview table of the x4 effect (and SE) from all msm_* objects you fit
# (assumes the objects exist in your workspace: msm_cl, msm_sl, msm_cn, msm_sn,
#  msm_clC, msm_slC, msm_cnC, msm_snC, msm_cnCI, msm_snCI)

effect_table <- data.frame(
  scenario = c(
    "1.1 constant_lin (no C)",
    "1.2 stepwise_lin (no C)",
    "1.3 constant_nonlin (no C)",
    "1.4 stepwise_nonlin (no C)",
    "2.1 constant_lin (with C1-3)",
    "2.2 stepwise_lin (with C1-3)",
    "2.3 constant_nonlin (with C1-3)",
    "2.4 stepwise_nonlin (with C1-3)",
    "3.1 constant_nonlin (with C + interactions)",
    "3.2 stepwise_nonlin (with C + interactions)"
  ),
  beta_x4 = c(
    coef(msm_cl)["x4"],
    coef(msm_sl)["x4"],
    coef(msm_cn)["x4"],
    coef(msm_sn)["x4"],
    coef(msm_clC)["x4"],
    coef(msm_slC)["x4"],
    coef(msm_cnC)["x4"],
    coef(msm_snC)["x4"],
    coef(msm_cnCI)["x4"],
    coef(msm_snCI)["x4"]
  ),
  se_x4 = c(
    summary(msm_cl)$coefficients["x4","Std. Error"],
    summary(msm_sl)$coefficients["x4","Std. Error"],
    summary(msm_cn)$coefficients["x4","Std. Error"],
    summary(msm_sn)$coefficients["x4","Std. Error"],
    summary(msm_clC)$coefficients["x4","Std. Error"],
    summary(msm_slC)$coefficients["x4","Std. Error"],
    summary(msm_cnC)$coefficients["x4","Std. Error"],
    summary(msm_snC)$coefficients["x4","Std. Error"],
    summary(msm_cnCI)$coefficients["x4","Std. Error"],
    summary(msm_snCI)$coefficients["x4","Std. Error"]
  )
)

effect_table$ci_low  <- effect_table$beta_x4 - 1.96 * effect_table$se_x4
effect_table$ci_high <- effect_table$beta_x4 + 1.96 * effect_table$se_x4

effect_table
