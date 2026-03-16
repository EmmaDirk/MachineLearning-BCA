# this script serves to test the function that residualise X and Y with respect to C
# step 1: define generate some testing data
# step 2: run the marginal regressions
# step 3: residualise the data
# step 4: run the marginal regressions on the residualised data
# step 5: check that the results are equal to the true parameters

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

# ------------------------------- 2) run regressions -------------------------------

lm1 <- lm(y2 ~ x1 + y1, data = data)
summary(lm1)

lm2 <- lm(x2 ~ x1 + y1, data = data)
summary(lm2)

# ------------------------------- 3) residualise the data -------------------------------

# A) linear residualiser
# residualise with only main effects
data_res_none <- residualise_panel_linearC(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 1
)

# residualise with all confounders
data_res_all <- residualise_panel_linearC(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 3
)

# residualise while ommiting c1
data_res_no_c1 <- residualise_panel_linearC(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = "c1",
  interaction_order = 3
)

# B) xgboost residualiser

# define a cheap tuning grid
xgb_tuning_grid <- expand.grid(
  eta = c(0.10, 0.20),
  max_depth = c(2, 3),
  min_child_weight = c(1, 5),
  subsample = 0.8,
  colsample_bytree = 0.8
)

# tune xgboost with all confounders
xgb_tuning_all <- tune_residualise_panel_xgb(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 3,
  tuning_grid = xgb_tuning_grid,
  cv_folds = 3,
  nrounds_max = 50,
  early_stopping_rounds = 5,
  nthread = 1,
  seed = 1234
)

# residualise using the tuned xgboost hyperparameters with all confounders
data_res_xgb_all <- residualise_panel_xgb(
  df = data,
  tuning = xgb_tuning_all,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 3,
  oof_folds = 2,
  nthread = 1,
  seed = 1234
)

# tune xgboost with no interaction confounders
# note: in the current xgboost residualiser, interaction_order is not used internally
# so this scenario is included to mirror the linear-model scenario naming
xgb_tuning_none <- tune_residualise_panel_xgb(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 1,
  tuning_grid = xgb_tuning_grid,
  cv_folds = 3,
  nrounds_max = 50,
  early_stopping_rounds = 5,
  nthread = 1,
  seed = 1234
)

# residualise using the tuned xgboost hyperparameters with no interaction confounders
# note: this will be equivalent to the all-confounder xgboost run unless the residualiser is changed
data_res_xgb_none <- residualise_panel_xgb(
  df = data,
  tuning = xgb_tuning_none,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = NULL,
  interaction_order = 1,
  oof_folds = 2,
  nthread = 1,
  seed = 1234
)

# tune xgboost while ommiting c1
xgb_tuning_no_c1 <- tune_residualise_panel_xgb(
  df = data,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = "c1",
  interaction_order = 3,
  tuning_grid = xgb_tuning_grid,
  cv_folds = 3,
  nrounds_max = 50,
  early_stopping_rounds = 5,
  nthread = 1,
  seed = 1234
)

# residualise using the tuned xgboost hyperparameters while ommiting c1
data_res_xgb_no_c1 <- residualise_panel_xgb(
  df = data,
  tuning = xgb_tuning_no_c1,
  k = 3,
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  exclude = "c1",
  interaction_order = 3,
  oof_folds = 2,
  nthread = 1,
  seed = 1234
)

# C) penalized regression method
# TBA) Penalized Regression (Lasso)

# TBA) TabPFN

# ------------------------------- 4) run regressions -------------------------------

# run regressions on residualised data with all confounders
# A) linear residualiser
lm3 <- lm(y2 ~ x1 + y1, data = data_res_all)
summary(lm3)

lm4 <- lm(x2 ~ x1 + y1, data = data_res_all)
summary(lm4)

# B) xgboost residualiser
xgblm1 <- lm(y2 ~ x1 + y1, data = data_res_xgb_all)
summary(xgblm1)

xgblm2 <- lm(x2 ~ x1 + y1, data = data_res_xgb_all)
summary(xgblm2)

# C) penalized regression method
# TBA) Penalized Regression (Lasso)

# TBA) TabPFN

# run residualised data with no interaction confounders
# A) linear residualiser
lm5 <- lm(y2 ~ x1 + y1, data = data_res_none)
summary(lm5)

lm6 <- lm(x2 ~ x1 + y1, data = data_res_none)
summary(lm6)

# B) xgboost residualiser
xgblm3 <- lm(y2 ~ x1 + y1, data = data_res_xgb_none)
summary(xgblm3)

xgblm4 <- lm(x2 ~ x1 + y1, data = data_res_xgb_none)
summary(xgblm4)

# C) penalized regression method
# TBA) Penalized Regression (Lasso)

# TBA) TabPFN

# run residualised data without c1
# A) linear residualiser
lm7 <- lm(y2 ~ x1 + y1, data = data_res_no_c1)
summary(lm7)

lm8 <- lm(x2 ~ x1 + y1, data = data_res_no_c1)
summary(lm8)

# B) xgboost residualiser
xgblm5 <- lm(y2 ~ x1 + y1, data = data_res_xgb_no_c1)
summary(xgblm5)

xgblm6 <- lm(x2 ~ x1 + y1, data = data_res_xgb_no_c1)
summary(xgblm6)

# TBA) lasso regression
# TBA) TabPFN

# ------------------------------- 5) compare to true parameters -------------------------------

# helper function to compute absolute deviations
extract_abs_dev <- function(model, method, true_x1, true_y1) {
  
  est <- coef(model)
  
  tibble(
    model = method,
    abs_dev_x1 = abs(unname(est["x1"]) - true_x1),
    abs_dev_y1 = abs(unname(est["y1"]) - true_y1)
  )
}

# collect results
deviation_results <- bind_rows(
  
  extract_abs_dev(lm3,    "linear_all",     Phi[2,1], Phi[2,2]),
  extract_abs_dev(lm5,    "linear_main",    Phi[2,1], Phi[2,2]),
  extract_abs_dev(lm7,    "linear_omit_c1", Phi[2,1], Phi[2,2]),
  extract_abs_dev(lm1,    "raw",            Phi[2,1], Phi[2,2]),
  
  extract_abs_dev(xgblm1, "xgb_all",        Phi[2,1], Phi[2,2]),
  extract_abs_dev(xgblm3, "xgb_main",       Phi[2,1], Phi[2,2]),
  extract_abs_dev(xgblm5, "xgb_omit_c1",    Phi[2,1], Phi[2,2])
)

# enforce model ordering
model_order <- c(
  "linear_all",
  "linear_main",
  "linear_omit_c1",
  "raw",
  "xgb_all",
  "xgb_main",
  "xgb_omit_c1"
)

deviation_results <- deviation_results %>%
  mutate(model = factor(model, levels = model_order)) %>%
  arrange(model) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  rename(
    Model = model,
    `Abs. deviation x1` = abs_dev_x1,
    `Abs. deviation y1` = abs_dev_y1
  )

# print table
deviation_results %>%
  kable(
    align = "lrr",
    caption = "Absolute deviations of estimated coefficients from the true parameters"
  )
