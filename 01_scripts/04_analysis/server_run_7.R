# this script is meant to analyze the data generated from the ferst department server run. 
# there are 24 scenarios:
#
# We always use:
# - T = 5
# - burn_in = 20
# - R^2_total = 0.05 
# - k = 5 with
# - omega11 = {1,  0,  0.05, 0.1,  0.15
#                , 1,  0.20, 0.25, 0.30
#                ,  ,     1, 0.35, 0.40
#                ,  ,      ,    1, 0.45}
# - Phi = {0.20, 0.00
#          0.10, 0.20}
# - Sigma = {1.00, 0.30
#            0.30, 1.00}
# 
# We test the following scenarios:
# 
# A) 10.000 observations
# 1) constant linear confounding effects
# 2) constant non-linear confounding effects, R^2_nl = 0.02
# 3) stepwise (rank order stable) non-linear confounding effects,
#    R^2_nl = 0.02, new_R^2 = 0.08, step_at = 3
# 4) stepwise (rank order unstable) non-linear confounding effects,
#    R^2_nl = 0.02, new_R^2 = 0.08, step_at = 3, lambda = 0.5
# 5) the same specification as 3, but we omit c1 and c2 
#    from every model except the true model.
# 6) the same specification as 4, but we omit c1 and c2
#    from every model except the true model.
#
# B) 1000 observations (1-6)
#
# C) 300 observations (1-6)
#
# D) 150 observations (1-6)
#
# For every scenario we fit:
# - the true model (CLPM with correct confounder order)
# - worst case model (CLPM with no adjustment)
# - the CLPM model with linear residualisation
# - the CLPM model with xgb residualisation
# - the CLPM model with enet residualisation
#
# - the RI-CLPM model without free loadings
# - the RI-CLPM model with free loadings
# - the RI-CLPM model with linear residualisation
# - the RI-CLPM model with xgb residualisation
# - the RI-CLPM model with enet residualisation
#
# - the DPM model without free loadings
# - the DPM model with free loadings
# - the DPM model with linear residualisation
# - the DPM model with xgb residualisation
# - the DPM model with enet residualisation
#
# In the non-linear scenarios we use both 2-way and 3-way interactions in the DGM.
# Therefore the correct confounder order for the true model is 3 in scenarios 2-6.
#
# We save one data frame per scenario x sample size combination.

# ============================================================
# libraries + setup
# ============================================================
library(here)
library(tidyverse)

# source plotting function
source(here("01_scripts", "04_analysis", "plotting.R"))

# ============================================================
# load data
# ============================================================

# N = 2000
dat_N2000_s1 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s01_N02000.rds"))
dat_N2000_s2 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s02_N02000.rds"))
dat_N2000_s3 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s03_N02000.rds"))
dat_N2000_s5 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s05_N02000.rds"))

# N = 1000
dat_N1000_s1 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s01_N01000.rds"))
dat_N1000_s2 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s02_N01000.rds"))
dat_N1000_s3 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s03_N01000.rds"))
dat_N1000_s5 <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s05_N01000.rds"))

# N = 300
dat_N300_s1  <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s01_N00300.rds"))
dat_N300_s2  <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s02_N00300.rds"))
dat_N300_s3  <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s03_N00300.rds"))
dat_N300_s5  <- readRDS(here("02_data", "02_server_runs", "07_run", "reduced_s05_N00300.rds"))

# ============================================================
# results
# ============================================================
out_s1_N2000 <- plot_overview_suite(dat_N2000_s1)

# ============================================================
# output scenario 1
# delta trajectory:
out_s1_N2000$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s1_N2000$plot_mse_x
out_s1_N2000$plot_mse_y
out_s1_N2000$plot_r2_x
out_s1_N2000$plot_r2_y

# convergence / proper solution
out_s1_N2000$plot_flag0
out_s1_N2000$plot_flag1
out_s1_N2000$plot_flag2
out_s1_N2000$plot_improper_reasons

# s.e. checks
out_s1_N2000$plot_se_ratio_clpm
out_s1_N2000$plot_se_ratio_riclpm
out_s1_N2000$plot_se_ratio_dpm
out_s1_N2000$plot_se_diff_clpm
out_s1_N2000$plot_se_diff_riclpm
out_s1_N2000$plot_se_diff_dpm

# clpm results
out_s1_N2000$plot_relbias_clpm
out_s1_N2000$plot_se_clpm
out_s1_N2000$plot_rmse_clpm
out_s1_N2000$plot_power_clpm

# ri-clpm results
out_s1_N2000$plot_relbias_riclpm
out_s1_N2000$plot_se_riclpm
out_s1_N2000$plot_rmse_riclpm
out_s1_N2000$plot_power_riclpm

# dpm results
out_s1_N2000$plot_relbias_dpm
out_s1_N2000$plot_se_dpm
out_s1_N2000$plot_rmse_dpm
out_s1_N2000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s1_N01000$plot_r2_x
out_s1_N01000$plot_r2_y
out_s1_N01000$plot_mse_x
out_s1_N01000$plot_mse_y

# convergence / proper solution
out_s1_N01000$plot_flag0
out_s1_N01000$plot_flag1
out_s1_N01000$plot_flag2
out_s1_N01000$plot_improper_reasons

# s.e. checks
out_s1_N01000$plot_se_ratio_clpm
out_s1_N01000$plot_se_ratio_riclpm
out_s1_N01000$plot_se_ratio_dpm
out_s1_N01000$plot_se_diff_clpm
out_s1_N01000$plot_se_diff_riclpm
out_s1_N01000$plot_se_diff_dpm

# clpm results
out_s1_N01000$plot_relbias_clpm
out_s1_N01000$plot_se_clpm
out_s1_N01000$plot_rmse_clpm
out_s1_N01000$plot_power_clpm

# ri-clpm results
out_s1_N01000$plot_relbias_riclpm
out_s1_N01000$plot_se_riclpm
out_s1_N01000$plot_rmse_riclpm
out_s1_N01000$plot_power_riclpm

# dpm results
out_s1_N01000$plot_relbias_dpm
out_s1_N01000$plot_se_dpm
out_s1_N01000$plot_rmse_dpm
out_s1_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s1_N00300$plot_r2_x
out_s1_N00300$plot_r2_y
out_s1_N00300$plot_mse_x
out_s1_N00300$plot_mse_y

# convergence / proper solution
out_s1_N00300$plot_flag0
out_s1_N00300$plot_flag1
out_s1_N00300$plot_flag2
out_s1_N00300$plot_improper_reasons

# s.e. checks
out_s1_N00300$plot_se_ratio_clpm
out_s1_N00300$plot_se_ratio_riclpm
out_s1_N00300$plot_se_ratio_dpm
out_s1_N00300$plot_se_diff_clpm
out_s1_N00300$plot_se_diff_riclpm
out_s1_N00300$plot_se_diff_dpm

# clpm results
out_s1_N00300$plot_relbias_clpm
out_s1_N00300$plot_se_clpm
out_s1_N00300$plot_rmse_clpm
out_s1_N00300$plot_power_clpm

# ri-clpm results
out_s1_N00300$plot_relbias_riclpm
out_s1_N00300$plot_se_riclpm
out_s1_N00300$plot_rmse_riclpm
out_s1_N00300$plot_power_riclpm

# dpm results
out_s1_N00300$plot_relbias_dpm
out_s1_N00300$plot_se_dpm
out_s1_N00300$plot_rmse_dpm
out_s1_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s1_N00150$plot_r2_x
out_s1_N00150$plot_r2_y
out_s1_N00150$plot_mse_x
out_s1_N00150$plot_mse_y

# convergence / proper solution
out_s1_N00150$plot_flag0
out_s1_N00150$plot_flag1
out_s1_N00150$plot_flag2
out_s1_N00150$plot_improper_reasons

# s.e. checks
out_s1_N00150$plot_se_ratio_clpm
out_s1_N00150$plot_se_ratio_riclpm
out_s1_N00150$plot_se_ratio_dpm
out_s1_N00150$plot_se_diff_clpm
out_s1_N00150$plot_se_diff_riclpm
out_s1_N00150$plot_se_diff_dpm

# clpm results
out_s1_N00150$plot_relbias_clpm
out_s1_N00150$plot_se_clpm
out_s1_N00150$plot_rmse_clpm
out_s1_N00150$plot_power_clpm

# ri-clpm results
out_s1_N00150$plot_relbias_riclpm
out_s1_N00150$plot_se_riclpm
out_s1_N00150$plot_rmse_riclpm
out_s1_N00150$plot_power_riclpm

# dpm results
out_s1_N00150$plot_relbias_dpm
out_s1_N00150$plot_se_dpm
out_s1_N00150$plot_rmse_dpm
out_s1_N00150$plot_power_dpm

# ============================================================
# output scenario 2
# delta trajectory:
out_s2_N00300$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s2_N10000$plot_r2_x
out_s2_N10000$plot_r2_y
out_s2_N10000$plot_mse_x
out_s2_N10000$plot_mse_y

# convergence / proper solution
out_s2_N10000$plot_flag0
out_s2_N10000$plot_flag1
out_s2_N10000$plot_flag2
out_s2_N10000$plot_improper_reasons

# s.e. checks
out_s2_N10000$plot_se_ratio_clpm
out_s2_N10000$plot_se_ratio_riclpm
out_s2_N10000$plot_se_ratio_dpm
out_s2_N10000$plot_se_diff_clpm
out_s2_N10000$plot_se_diff_riclpm
out_s2_N10000$plot_se_diff_dpm

# clpm results
out_s2_N10000$plot_relbias_clpm
out_s2_N10000$plot_se_clpm
out_s2_N10000$plot_rmse_clpm
out_s2_N10000$plot_power_clpm

# ri-clpm results
out_s2_N10000$plot_relbias_riclpm
out_s2_N10000$plot_se_riclpm
out_s2_N10000$plot_rmse_riclpm
out_s2_N10000$plot_power_riclpm

# dpm results
out_s2_N10000$plot_relbias_dpm
out_s2_N10000$plot_se_dpm
out_s2_N10000$plot_rmse_dpm
out_s2_N10000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s2_N01000$plot_r2_x
out_s2_N01000$plot_r2_y
out_s2_N01000$plot_mse_x
out_s2_N01000$plot_mse_y

# convergence / proper solution
out_s2_N01000$plot_flag0
out_s2_N01000$plot_flag1
out_s2_N01000$plot_flag2
out_s2_N01000$plot_improper_reasons

# s.e. checks
out_s2_N01000$plot_se_ratio_clpm
out_s2_N01000$plot_se_ratio_riclpm
out_s2_N01000$plot_se_ratio_dpm
out_s2_N01000$plot_se_diff_clpm
out_s2_N01000$plot_se_diff_riclpm
out_s2_N01000$plot_se_diff_dpm

# clpm results
out_s2_N01000$plot_relbias_clpm
out_s2_N01000$plot_se_clpm
out_s2_N01000$plot_rmse_clpm
out_s2_N01000$plot_power_clpm

# ri-clpm results
out_s2_N01000$plot_relbias_riclpm
out_s2_N01000$plot_se_riclpm
out_s2_N01000$plot_rmse_riclpm
out_s2_N01000$plot_power_riclpm

# dpm results
out_s2_N01000$plot_relbias_dpm
out_s2_N01000$plot_se_dpm
out_s2_N01000$plot_rmse_dpm
out_s2_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s2_N00300$plot_r2_x
out_s2_N00300$plot_r2_y
out_s2_N00300$plot_mse_x
out_s2_N00300$plot_mse_y

# convergence / proper solution
out_s2_N00300$plot_flag0
out_s2_N00300$plot_flag1
out_s2_N00300$plot_flag2
out_s2_N00300$plot_improper_reasons

# s.e. checks
out_s2_N00300$plot_se_ratio_clpm
out_s2_N00300$plot_se_ratio_riclpm
out_s2_N00300$plot_se_ratio_dpm
out_s2_N00300$plot_se_diff_clpm
out_s2_N00300$plot_se_diff_riclpm
out_s2_N00300$plot_se_diff_dpm

# clpm results
out_s2_N00300$plot_relbias_clpm
out_s2_N00300$plot_se_clpm
out_s2_N00300$plot_rmse_clpm
out_s2_N00300$plot_power_clpm

# ri-clpm results
out_s2_N00300$plot_relbias_riclpm
out_s2_N00300$plot_se_riclpm
out_s2_N00300$plot_rmse_riclpm
out_s2_N00300$plot_power_riclpm

# dpm results
out_s2_N00300$plot_relbias_dpm
out_s2_N00300$plot_se_dpm
out_s2_N00300$plot_rmse_dpm
out_s2_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s2_N00150$plot_r2_x
out_s2_N00150$plot_r2_y
out_s2_N00150$plot_mse_x
out_s2_N00150$plot_mse_y

# convergence / proper solution
out_s2_N00150$plot_flag0
out_s2_N00150$plot_flag1
out_s2_N00150$plot_flag2
out_s2_N00150$plot_improper_reasons

# s.e. checks
out_s2_N00150$plot_se_ratio_clpm
out_s2_N00150$plot_se_ratio_riclpm
out_s2_N00150$plot_se_ratio_dpm
out_s2_N00150$plot_se_diff_clpm
out_s2_N00150$plot_se_diff_riclpm
out_s2_N00150$plot_se_diff_dpm

# clpm results
out_s2_N00150$plot_relbias_clpm
out_s2_N00150$plot_se_clpm
out_s2_N00150$plot_rmse_clpm
out_s2_N00150$plot_power_clpm

# ri-clpm results
out_s2_N00150$plot_relbias_riclpm
out_s2_N00150$plot_se_riclpm
out_s2_N00150$plot_rmse_riclpm
out_s2_N00150$plot_power_riclpm

# dpm results
out_s2_N00150$plot_relbias_dpm
out_s2_N00150$plot_se_dpm
out_s2_N00150$plot_rmse_dpm
out_s2_N00150$plot_power_dpm

# ============================================================
# output scenario 3
# delta trajectory:
out_s3_N00300$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s3_N10000$plot_r2_x
out_s3_N10000$plot_r2_y
out_s3_N10000$plot_mse_x
out_s3_N10000$plot_mse_y

# convergence / proper solution
out_s3_N10000$plot_flag0
out_s3_N10000$plot_flag1
out_s3_N10000$plot_flag2
out_s3_N10000$plot_improper_reasons

# s.e. checks
out_s3_N10000$plot_se_ratio_clpm
out_s3_N10000$plot_se_ratio_riclpm
out_s3_N10000$plot_se_ratio_dpm
out_s3_N10000$plot_se_diff_clpm
out_s3_N10000$plot_se_diff_riclpm
out_s3_N10000$plot_se_diff_dpm

# clpm results
out_s3_N10000$plot_relbias_clpm
out_s3_N10000$plot_se_clpm
out_s3_N10000$plot_rmse_clpm
out_s3_N10000$plot_power_clpm

# ri-clpm results
out_s3_N10000$plot_relbias_riclpm
out_s3_N10000$plot_se_riclpm
out_s3_N10000$plot_rmse_riclpm
out_s3_N10000$plot_power_riclpm

# dpm results
out_s3_N10000$plot_relbias_dpm
out_s3_N10000$plot_se_dpm
out_s3_N10000$plot_rmse_dpm
out_s3_N10000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s3_N01000$plot_r2_x
out_s3_N01000$plot_r2_y
out_s3_N01000$plot_mse_x
out_s3_N01000$plot_mse_y

# convergence / proper solution
out_s3_N01000$plot_flag0
out_s3_N01000$plot_flag1
out_s3_N01000$plot_flag2
out_s3_N01000$plot_improper_reasons

# s.e. checks
out_s3_N01000$plot_se_ratio_clpm
out_s3_N01000$plot_se_ratio_riclpm
out_s3_N01000$plot_se_ratio_dpm
out_s3_N01000$plot_se_diff_clpm
out_s3_N01000$plot_se_diff_riclpm
out_s3_N01000$plot_se_diff_dpm

# clpm results
out_s3_N01000$plot_relbias_clpm
out_s3_N01000$plot_se_clpm
out_s3_N01000$plot_rmse_clpm
out_s3_N01000$plot_power_clpm

# ri-clpm results
out_s3_N01000$plot_relbias_riclpm
out_s3_N01000$plot_se_riclpm
out_s3_N01000$plot_rmse_riclpm
out_s3_N01000$plot_power_riclpm

# dpm results
out_s3_N01000$plot_relbias_dpm
out_s3_N01000$plot_se_dpm
out_s3_N01000$plot_rmse_dpm
out_s3_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s3_N00300$plot_r2_x
out_s3_N00300$plot_r2_y
out_s3_N00300$plot_mse_x
out_s3_N00300$plot_mse_y

# convergence / proper solution
out_s3_N00300$plot_flag0
out_s3_N00300$plot_flag1
out_s3_N00300$plot_flag2
out_s3_N00300$plot_improper_reasons

# s.e. checks
out_s3_N00300$plot_se_ratio_clpm
out_s3_N00300$plot_se_ratio_riclpm
out_s3_N00300$plot_se_ratio_dpm
out_s3_N00300$plot_se_diff_clpm
out_s3_N00300$plot_se_diff_riclpm
out_s3_N00300$plot_se_diff_dpm

# clpm results
out_s3_N00300$plot_relbias_clpm
out_s3_N00300$plot_se_clpm
out_s3_N00300$plot_rmse_clpm
out_s3_N00300$plot_power_clpm

# ri-clpm results
out_s3_N00300$plot_relbias_riclpm
out_s3_N00300$plot_se_riclpm
out_s3_N00300$plot_rmse_riclpm
out_s3_N00300$plot_power_riclpm

# dpm results
out_s3_N00300$plot_relbias_dpm
out_s3_N00300$plot_se_dpm
out_s3_N00300$plot_rmse_dpm
out_s3_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s3_N00150$plot_r2_x
out_s3_N00150$plot_r2_y
out_s3_N00150$plot_mse_x
out_s3_N00150$plot_mse_y

# convergence / proper solution
out_s3_N00150$plot_flag0
out_s3_N00150$plot_flag1
out_s3_N00150$plot_flag2
out_s3_N00150$plot_improper_reasons

# s.e. checks
out_s3_N00150$plot_se_ratio_clpm
out_s3_N00150$plot_se_ratio_riclpm
out_s3_N00150$plot_se_ratio_dpm
out_s3_N00150$plot_se_diff_clpm
out_s3_N00150$plot_se_diff_riclpm
out_s3_N00150$plot_se_diff_dpm

# clpm results
out_s3_N00150$plot_relbias_clpm
out_s3_N00150$plot_se_clpm
out_s3_N00150$plot_rmse_clpm
out_s3_N00150$plot_power_clpm

# ri-clpm results
out_s3_N00150$plot_relbias_riclpm
out_s3_N00150$plot_se_riclpm
out_s3_N00150$plot_rmse_riclpm
out_s3_N00150$plot_power_riclpm

# dpm results
out_s3_N00150$plot_relbias_dpm
out_s3_N00150$plot_se_dpm
out_s3_N00150$plot_rmse_dpm
out_s3_N00150$plot_power_dpm

# ============================================================
# output scenario 4
# delta trajectory:
out_s4_N00300$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s4_N10000$plot_r2_x
out_s4_N10000$plot_r2_y
out_s4_N10000$plot_mse_x
out_s4_N10000$plot_mse_y

# convergence / proper solution
out_s4_N10000$plot_flag0
out_s4_N10000$plot_flag1
out_s4_N10000$plot_flag2
out_s4_N10000$plot_improper_reasons

# s.e. checks
out_s4_N10000$plot_se_ratio_clpm
out_s4_N10000$plot_se_ratio_riclpm
out_s4_N10000$plot_se_ratio_dpm
out_s4_N10000$plot_se_diff_clpm
out_s4_N10000$plot_se_diff_riclpm
out_s4_N10000$plot_se_diff_dpm

# clpm results
out_s4_N10000$plot_relbias_clpm
out_s4_N10000$plot_se_clpm
out_s4_N10000$plot_rmse_clpm
out_s4_N10000$plot_power_clpm

# ri-clpm results
out_s4_N10000$plot_relbias_riclpm
out_s4_N10000$plot_se_riclpm
out_s4_N10000$plot_rmse_riclpm
out_s4_N10000$plot_power_riclpm

# dpm results
out_s4_N10000$plot_relbias_dpm
out_s4_N10000$plot_se_dpm
out_s4_N10000$plot_rmse_dpm
out_s4_N10000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s4_N01000$plot_r2_x
out_s4_N01000$plot_r2_y
out_s4_N01000$plot_mse_x
out_s4_N01000$plot_mse_y

# convergence / proper solution
out_s4_N01000$plot_flag0
out_s4_N01000$plot_flag1
out_s4_N01000$plot_flag2
out_s4_N01000$plot_improper_reasons

# s.e. checks
out_s4_N01000$plot_se_ratio_clpm
out_s4_N01000$plot_se_ratio_riclpm
out_s4_N01000$plot_se_ratio_dpm
out_s4_N01000$plot_se_diff_clpm
out_s4_N01000$plot_se_diff_riclpm
out_s4_N01000$plot_se_diff_dpm

# clpm results
out_s4_N01000$plot_relbias_clpm
out_s4_N01000$plot_se_clpm
out_s4_N01000$plot_rmse_clpm
out_s4_N01000$plot_power_clpm

# ri-clpm results
out_s4_N01000$plot_relbias_riclpm
out_s4_N01000$plot_se_riclpm
out_s4_N01000$plot_rmse_riclpm
out_s4_N01000$plot_power_riclpm

# dpm results
out_s4_N01000$plot_relbias_dpm
out_s4_N01000$plot_se_dpm
out_s4_N01000$plot_rmse_dpm
out_s4_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s4_N00300$plot_r2_x
out_s4_N00300$plot_r2_y
out_s4_N00300$plot_mse_x
out_s4_N00300$plot_mse_y

# convergence / proper solution
out_s4_N00300$plot_flag0
out_s4_N00300$plot_flag1
out_s4_N00300$plot_flag2
out_s4_N00300$plot_improper_reasons

# s.e. checks
out_s4_N00300$plot_se_ratio_clpm
out_s4_N00300$plot_se_ratio_riclpm
out_s4_N00300$plot_se_ratio_dpm
out_s4_N00300$plot_se_diff_clpm
out_s4_N00300$plot_se_diff_riclpm
out_s4_N00300$plot_se_diff_dpm

# clpm results
out_s4_N00300$plot_relbias_clpm
out_s4_N00300$plot_se_clpm
out_s4_N00300$plot_rmse_clpm
out_s4_N00300$plot_power_clpm

# ri-clpm results
out_s4_N00300$plot_relbias_riclpm
out_s4_N00300$plot_se_riclpm
out_s4_N00300$plot_rmse_riclpm
out_s4_N00300$plot_power_riclpm

# dpm results
out_s4_N00300$plot_relbias_dpm
out_s4_N00300$plot_se_dpm
out_s4_N00300$plot_rmse_dpm
out_s4_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s4_N00150$plot_r2_x
out_s4_N00150$plot_r2_y
out_s4_N00150$plot_mse_x
out_s4_N00150$plot_mse_y

# convergence / proper solution
out_s4_N00150$plot_flag0
out_s4_N00150$plot_flag1
out_s4_N00150$plot_flag2
out_s4_N00150$plot_improper_reasons

# s.e. checks
out_s4_N00150$plot_se_ratio_clpm
out_s4_N00150$plot_se_ratio_riclpm
out_s4_N00150$plot_se_ratio_dpm
out_s4_N00150$plot_se_diff_clpm
out_s4_N00150$plot_se_diff_riclpm
out_s4_N00150$plot_se_diff_dpm

# clpm results
out_s4_N00150$plot_relbias_clpm
out_s4_N00150$plot_se_clpm
out_s4_N00150$plot_rmse_clpm
out_s4_N00150$plot_power_clpm

# ri-clpm results
out_s4_N00150$plot_relbias_riclpm
out_s4_N00150$plot_se_riclpm
out_s4_N00150$plot_rmse_riclpm
out_s4_N00150$plot_power_riclpm

# dpm results
out_s4_N00150$plot_relbias_dpm
out_s4_N00150$plot_se_dpm
out_s4_N00150$plot_rmse_dpm
out_s4_N00150$plot_power_dpm

# ============================================================
# output scenario 5
# delta trajectory:
out_s5_N00300$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s5_N10000$plot_r2_x
out_s5_N10000$plot_r2_y
out_s5_N10000$plot_mse_x
out_s5_N10000$plot_mse_y

# convergence / proper solution
out_s5_N10000$plot_flag0
out_s5_N10000$plot_flag1
out_s5_N10000$plot_flag2
out_s5_N10000$plot_improper_reasons

# s.e. checks
out_s5_N10000$plot_se_ratio_clpm
out_s5_N10000$plot_se_ratio_riclpm
out_s5_N10000$plot_se_ratio_dpm
out_s5_N10000$plot_se_diff_clpm
out_s5_N10000$plot_se_diff_riclpm
out_s5_N10000$plot_se_diff_dpm

# clpm results
out_s5_N10000$plot_relbias_clpm
out_s5_N10000$plot_se_clpm
out_s5_N10000$plot_rmse_clpm
out_s5_N10000$plot_power_clpm

# ri-clpm results
out_s5_N10000$plot_relbias_riclpm
out_s5_N10000$plot_se_riclpm
out_s5_N10000$plot_rmse_riclpm
out_s5_N10000$plot_power_riclpm

# dpm results
out_s5_N10000$plot_relbias_dpm
out_s5_N10000$plot_se_dpm
out_s5_N10000$plot_rmse_dpm
out_s5_N10000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s5_N01000$plot_r2_x
out_s5_N01000$plot_r2_y
out_s5_N01000$plot_mse_x
out_s5_N01000$plot_mse_y

# convergence / proper solution
out_s5_N01000$plot_flag0
out_s5_N01000$plot_flag1
out_s5_N01000$plot_flag2
out_s5_N01000$plot_improper_reasons

# s.e. checks
out_s5_N01000$plot_se_ratio_clpm
out_s5_N01000$plot_se_ratio_riclpm
out_s5_N01000$plot_se_ratio_dpm
out_s5_N01000$plot_se_diff_clpm
out_s5_N01000$plot_se_diff_riclpm
out_s5_N01000$plot_se_diff_dpm

# clpm results
out_s5_N01000$plot_relbias_clpm
out_s5_N01000$plot_se_clpm
out_s5_N01000$plot_rmse_clpm
out_s5_N01000$plot_power_clpm

# ri-clpm results
out_s5_N01000$plot_relbias_riclpm
out_s5_N01000$plot_se_riclpm
out_s5_N01000$plot_rmse_riclpm
out_s5_N01000$plot_power_riclpm

# dpm results
out_s5_N01000$plot_relbias_dpm
out_s5_N01000$plot_se_dpm
out_s5_N01000$plot_rmse_dpm
out_s5_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s5_N00300$plot_r2_x
out_s5_N00300$plot_r2_y
out_s5_N00300$plot_mse_x
out_s5_N00300$plot_mse_y

# convergence / proper solution
out_s5_N00300$plot_flag0
out_s5_N00300$plot_flag1
out_s5_N00300$plot_flag2
out_s5_N00300$plot_improper_reasons

# s.e. checks
out_s5_N00300$plot_se_ratio_clpm
out_s5_N00300$plot_se_ratio_riclpm
out_s5_N00300$plot_se_ratio_dpm
out_s5_N00300$plot_se_diff_clpm
out_s5_N00300$plot_se_diff_riclpm
out_s5_N00300$plot_se_diff_dpm

# clpm results
out_s5_N00300$plot_relbias_clpm
out_s5_N00300$plot_se_clpm
out_s5_N00300$plot_rmse_clpm
out_s5_N00300$plot_power_clpm

# ri-clpm results
out_s5_N00300$plot_relbias_riclpm
out_s5_N00300$plot_se_riclpm
out_s5_N00300$plot_rmse_riclpm
out_s5_N00300$plot_power_riclpm

# dpm results
out_s5_N00300$plot_relbias_dpm
out_s5_N00300$plot_se_dpm
out_s5_N00300$plot_rmse_dpm
out_s5_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s5_N00150$plot_r2_x
out_s5_N00150$plot_r2_y
out_s5_N00150$plot_mse_x
out_s5_N00150$plot_mse_y

# convergence / proper solution
out_s5_N00150$plot_flag0
out_s5_N00150$plot_flag1
out_s5_N00150$plot_flag2
out_s5_N00150$plot_improper_reasons

# s.e. checks
out_s5_N00150$plot_se_ratio_clpm
out_s5_N00150$plot_se_ratio_riclpm
out_s5_N00150$plot_se_ratio_dpm
out_s5_N00150$plot_se_diff_clpm
out_s5_N00150$plot_se_diff_riclpm
out_s5_N00150$plot_se_diff_dpm

# clpm results
out_s5_N00150$plot_relbias_clpm
out_s5_N00150$plot_se_clpm
out_s5_N00150$plot_rmse_clpm
out_s5_N00150$plot_power_clpm

# ri-clpm results
out_s5_N00150$plot_relbias_riclpm
out_s5_N00150$plot_se_riclpm
out_s5_N00150$plot_rmse_riclpm
out_s5_N00150$plot_power_riclpm

# dpm results
out_s5_N00150$plot_relbias_dpm
out_s5_N00150$plot_se_dpm
out_s5_N00150$plot_rmse_dpm
out_s5_N00150$plot_power_dpm

# ============================================================
# output scenario 6
# delta trajectory:
out_s6_N00300$plot_true_delta
# ============================================================

# ------------------------ N = 10000 -------------------------

# ml diagnostics
out_s6_N10000$plot_r2_x
out_s6_N10000$plot_r2_y
out_s6_N10000$plot_mse_x
out_s6_N10000$plot_mse_y

# convergence / proper solution
out_s6_N10000$plot_flag0
out_s6_N10000$plot_flag1
out_s6_N10000$plot_flag2
out_s6_N10000$plot_improper_reasons

# s.e. checks
out_s6_N10000$plot_se_ratio_clpm
out_s6_N10000$plot_se_ratio_riclpm
out_s6_N10000$plot_se_ratio_dpm
out_s6_N10000$plot_se_diff_clpm
out_s6_N10000$plot_se_diff_riclpm
out_s6_N10000$plot_se_diff_dpm

# clpm results
out_s6_N10000$plot_relbias_clpm
out_s6_N10000$plot_se_clpm
out_s6_N10000$plot_rmse_clpm
out_s6_N10000$plot_power_clpm

# ri-clpm results
out_s6_N10000$plot_relbias_riclpm
out_s6_N10000$plot_se_riclpm
out_s6_N10000$plot_rmse_riclpm
out_s6_N10000$plot_power_riclpm

# dpm results
out_s6_N10000$plot_relbias_dpm
out_s6_N10000$plot_se_dpm
out_s6_N10000$plot_rmse_dpm
out_s6_N10000$plot_power_dpm

# ------------------------ N = 1000 --------------------------

# ml diagnostics
out_s6_N01000$plot_r2_x
out_s6_N01000$plot_r2_y
out_s6_N01000$plot_mse_x
out_s6_N01000$plot_mse_y

# convergence / proper solution
out_s6_N01000$plot_flag0
out_s6_N01000$plot_flag1
out_s6_N01000$plot_flag2
out_s6_N01000$plot_improper_reasons

# s.e. checks
out_s6_N01000$plot_se_ratio_clpm
out_s6_N01000$plot_se_ratio_riclpm
out_s6_N01000$plot_se_ratio_dpm
out_s6_N01000$plot_se_diff_clpm
out_s6_N01000$plot_se_diff_riclpm
out_s6_N01000$plot_se_diff_dpm

# clpm results
out_s6_N01000$plot_relbias_clpm
out_s6_N01000$plot_se_clpm
out_s6_N01000$plot_rmse_clpm
out_s6_N01000$plot_power_clpm

# ri-clpm results
out_s6_N01000$plot_relbias_riclpm
out_s6_N01000$plot_se_riclpm
out_s6_N01000$plot_rmse_riclpm
out_s6_N01000$plot_power_riclpm

# dpm results
out_s6_N01000$plot_relbias_dpm
out_s6_N01000$plot_se_dpm
out_s6_N01000$plot_rmse_dpm
out_s6_N01000$plot_power_dpm

# ------------------------ N = 300 ---------------------------

# ml diagnostics
out_s6_N00300$plot_r2_x
out_s6_N00300$plot_r2_y
out_s6_N00300$plot_mse_x
out_s6_N00300$plot_mse_y

# convergence / proper solution
out_s6_N00300$plot_flag0
out_s6_N00300$plot_flag1
out_s6_N00300$plot_flag2
out_s6_N00300$plot_improper_reasons

# s.e. checks
out_s6_N00300$plot_se_ratio_clpm
out_s6_N00300$plot_se_ratio_riclpm
out_s6_N00300$plot_se_ratio_dpm
out_s6_N00300$plot_se_diff_clpm
out_s6_N00300$plot_se_diff_riclpm
out_s6_N00300$plot_se_diff_dpm

# clpm results
out_s6_N00300$plot_relbias_clpm
out_s6_N00300$plot_se_clpm
out_s6_N00300$plot_rmse_clpm
out_s6_N00300$plot_power_clpm

# ri-clpm results
out_s6_N00300$plot_relbias_riclpm
out_s6_N00300$plot_se_riclpm
out_s6_N00300$plot_rmse_riclpm
out_s6_N00300$plot_power_riclpm

# dpm results
out_s6_N00300$plot_relbias_dpm
out_s6_N00300$plot_se_dpm
out_s6_N00300$plot_rmse_dpm
out_s6_N00300$plot_power_dpm

# ------------------------ N = 150 ---------------------------

# ml diagnostics
out_s6_N00150$plot_r2_x
out_s6_N00150$plot_r2_y
out_s6_N00150$plot_mse_x
out_s6_N00150$plot_mse_y

# convergence / proper solution
out_s6_N00150$plot_flag0
out_s6_N00150$plot_flag1
out_s6_N00150$plot_flag2
out_s6_N00150$plot_improper_reasons

# s.e. checks
out_s6_N00150$plot_se_ratio_clpm
out_s6_N00150$plot_se_ratio_riclpm
out_s6_N00150$plot_se_ratio_dpm
out_s6_N00150$plot_se_diff_clpm
out_s6_N00150$plot_se_diff_riclpm
out_s6_N00150$plot_se_diff_dpm

# clpm results
out_s6_N00150$plot_relbias_clpm
out_s6_N00150$plot_se_clpm
out_s6_N00150$plot_rmse_clpm
out_s6_N00150$plot_power_clpm

# ri-clpm results
out_s6_N00150$plot_relbias_riclpm
out_s6_N00150$plot_se_riclpm
out_s6_N00150$plot_rmse_riclpm
out_s6_N00150$plot_power_riclpm

# dpm results
out_s6_N00150$plot_relbias_dpm
out_s6_N00150$plot_se_dpm
out_s6_N00150$plot_rmse_dpm
out_s6_N00150$plot_power_dpm