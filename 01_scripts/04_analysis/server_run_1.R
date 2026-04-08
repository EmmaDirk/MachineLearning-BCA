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

source(here("01_scripts", "02_development", "06_thesis_v3", "00_packages.R"))

# source plotting function
source(here("01_scripts", "02_development", "06_thesis_v3", "12_plotting.R"))

# ============================================================
# load data
# ============================================================

# N = 10000
dat_N10000_s1 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s01_N10000.rds"))
dat_N10000_s2 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s02_N10000.rds"))
dat_N10000_s3 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s03_N10000.rds"))
dat_N10000_s4 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s04_N10000.rds"))
dat_N10000_s5 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s05_N10000.rds"))
dat_N10000_s6 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s06_N10000.rds"))

# N = 1000
dat_N01000_s1 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s01_N01000.rds"))
dat_N01000_s2 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s02_N01000.rds"))
dat_N01000_s3 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s03_N01000.rds"))
dat_N01000_s4 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s04_N01000.rds"))
dat_N01000_s5 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s05_N01000.rds"))
dat_N01000_s6 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s06_N01000.rds"))

# N = 300
dat_N00300_s1 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s01_N00300.rds"))
dat_N00300_s2 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s02_N00300.rds"))
dat_N00300_s3 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s03_N00300.rds"))
dat_N00300_s4 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s04_N00300.rds"))
dat_N00300_s5 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s05_N00300.rds"))
dat_N00300_s6 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s06_N00300.rds"))

# N = 150
dat_N00150_s1 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s01_N00150.rds"))
dat_N00150_s2 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s02_N00150.rds"))
dat_N00150_s3 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s03_N00150.rds"))
dat_N00150_s4 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s04_N00150.rds"))
dat_N00150_s5 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s05_N00150.rds"))
dat_N00150_s6 <- readRDS(here("02_data", "02_server_runs", "01_run", "test_s06_N00150.rds"))

# ============================================================
# results
# ============================================================
out_s1_N10000 <- plot_engine_results(dat_N10000_s1)
out_s1_N01000 <- plot_engine_results(dat_N01000_s1)
out_s1_N00300 <- plot_engine_results(dat_N00300_s1)
out_s1_N00150 <- plot_engine_results(dat_N00150_s1)

out_s2_N10000 <- plot_engine_results(dat_N10000_s2)
out_s2_N01000 <- plot_engine_results(dat_N01000_s2)
out_s2_N00300 <- plot_engine_results(dat_N00300_s2)
out_s2_N00150 <- plot_engine_results(dat_N00150_s2)

out_s3_N10000 <- plot_engine_results(dat_N10000_s3)
out_s3_N01000 <- plot_engine_results(dat_N01000_s3)
out_s3_N00300 <- plot_engine_results(dat_N00300_s3)
out_s3_N00150 <- plot_engine_results(dat_N00150_s3)

out_s4_N10000 <- plot_engine_results(dat_N10000_s4)
out_s4_N01000 <- plot_engine_results(dat_N01000_s4)
out_s4_N00300 <- plot_engine_results(dat_N00300_s4)
out_s4_N00150 <- plot_engine_results(dat_N00150_s4)

out_s5_N10000 <- plot_engine_results(dat_N10000_s5)
out_s5_N01000 <- plot_engine_results(dat_N01000_s5)
out_s5_N00300 <- plot_engine_results(dat_N00300_s5)
out_s5_N00150 <- plot_engine_results(dat_N00150_s5)

out_s6_N10000 <- plot_engine_results(dat_N10000_s6)
out_s6_N01000 <- plot_engine_results(dat_N01000_s6)
out_s6_N00300 <- plot_engine_results(dat_N00300_s6)
out_s6_N00150 <- plot_engine_results(dat_N00150_s6)

# ============================================================
# Scenario 1
# constant linear confounding effects
# ============================================================

# ---------------------- 10000 -------------------------------

# results CLPM
out_s1_N10000$plot_relbias_clpm
out_s1_N10000$plot_se_clpm
out_s1_N10000$plot_rmse_clpm
out_s1_N10000$plot_power_clpm

# results RI-CLPM
out_s1_N10000$plot_relbias_riclpm
out_s1_N10000$plot_se_riclpm
out_s1_N10000$plot_rmse_riclpm
out_s1_N10000$plot_power_riclpm

# results DPM
out_s1_N10000$plot_relbias_dpm
out_s1_N10000$plot_se_dpm
out_s1_N10000$plot_rmse_dpm
out_s1_N10000$plot_power_dpm

# convergence checks
print(out_s1_N10000$flag0_df, n = 100)
print(out_s1_N10000$flag1_df, n = 100)
print(out_s1_N10000$flag2_df, n = 100)

# standard error checks
out_s1_N10000$plot_se_ratio
out_s1_N10000$plot_se_diff

# ML diagnostics
out_s1_N10000$plot_r2
out_s1_N10000$plot_mse

# ---------------------- 01000 -------------------------------

# results CLPM
out_s1_N01000$plot_relbias_clpm
out_s1_N01000$plot_se_clpm
out_s1_N01000$plot_rmse_clpm
out_s1_N01000$plot_power_clpm

# results RI-CLPM
out_s1_N01000$plot_relbias_riclpm
out_s1_N01000$plot_se_riclpm
out_s1_N01000$plot_rmse_riclpm
out_s1_N01000$plot_power_riclpm

# results DPM
out_s1_N01000$plot_relbias_dpm
out_s1_N01000$plot_se_dpm
out_s1_N01000$plot_rmse_dpm
out_s1_N01000$plot_power_dpm

# convergence checks
print(out_s1_N01000$flag0_df, n = 100)
print(out_s1_N01000$flag1_df, n = 100)
print(out_s1_N01000$flag2_df, n = 100)

# standard error checks
out_s1_N01000$plot_se_ratio
out_s1_N01000$plot_se_diff

# ML diagnostics
out_s1_N01000$plot_r2
out_s1_N01000$plot_mse

# ---------------------- 00300 -------------------------------

# results CLPM
out_s1_N00300$plot_relbias_clpm
out_s1_N00300$plot_se_clpm
out_s1_N00300$plot_rmse_clpm
out_s1_N00300$plot_power_clpm

# results RI-CLPM
out_s1_N00300$plot_relbias_riclpm
out_s1_N00300$plot_se_riclpm
out_s1_N00300$plot_rmse_riclpm
out_s1_N00300$plot_power_riclpm

# results DPM
out_s1_N00300$plot_relbias_dpm
out_s1_N00300$plot_se_dpm
out_s1_N00300$plot_rmse_dpm
out_s1_N00300$plot_power_dpm

# convergence checks
print(out_s1_N00300$flag0_df, n = 100)
print(out_s1_N00300$flag1_df, n = 100)
print(out_s1_N00300$flag2_df, n = 100)

# standard error checks
out_s1_N00300$plot_se_ratio
out_s1_N00300$plot_se_diff

# ML diagnostics
out_s1_N00300$plot_r2
out_s1_N00300$plot_mse

# ---------------------- 00150 -------------------------------

# results CLPM
out_s1_N00150$plot_relbias_clpm
out_s1_N00150$plot_se_clpm
out_s1_N00150$plot_rmse_clpm
out_s1_N00150$plot_power_clpm

# results RI-CLPM
out_s1_N00150$plot_relbias_riclpm
out_s1_N00150$plot_se_riclpm
out_s1_N00150$plot_rmse_riclpm
out_s1_N00150$plot_power_riclpm

# results DPM
out_s1_N00150$plot_relbias_dpm
out_s1_N00150$plot_se_dpm
out_s1_N00150$plot_rmse_dpm
out_s1_N00150$plot_power_dpm

# convergence checks
print(out_s1_N00150$flag0_df, n = 100)
print(out_s1_N00150$flag1_df, n = 100)
print(out_s1_N00150$flag2_df, n = 100)

# standard error checks
out_s1_N00150$plot_se_ratio
out_s1_N00150$plot_se_diff

# ML diagnostics
out_s1_N00150$plot_r2
out_s1_N00150$plot_mse


# ============================================================
# Scenario 2
# constant non-linear confounding effects
# ============================================================

# ---------------------- 10000 -------------------------------

# results CLPM
out_s2_N10000$plot_relbias_clpm
out_s2_N10000$plot_se_clpm
out_s2_N10000$plot_rmse_clpm
out_s2_N10000$plot_power_clpm

# results RI-CLPM
out_s2_N10000$plot_relbias_riclpm
out_s2_N10000$plot_se_riclpm
out_s2_N10000$plot_rmse_riclpm
out_s2_N10000$plot_power_riclpm

# results DPM
out_s2_N10000$plot_relbias_dpm
out_s2_N10000$plot_se_dpm
out_s2_N10000$plot_rmse_dpm
out_s2_N10000$plot_power_dpm

# convergence checks
print(out_s2_N10000$flag0_df, n = 100)
print(out_s2_N10000$flag1_df, n = 100)
print(out_s2_N10000$flag2_df, n = 100)

# standard error checks
out_s2_N10000$plot_se_ratio
out_s2_N10000$plot_se_diff

# ML diagnostics
out_s2_N10000$plot_r2
out_s2_N10000$plot_mse


# ---------------------- 01000 -------------------------------

# results CLPM
out_s2_N01000$plot_relbias_clpm
out_s2_N01000$plot_se_clpm
out_s2_N01000$plot_rmse_clpm
out_s2_N01000$plot_power_clpm

# results RI-CLPM
out_s2_N01000$plot_relbias_riclpm
out_s2_N01000$plot_se_riclpm
out_s2_N01000$plot_rmse_riclpm
out_s2_N01000$plot_power_riclpm

# results DPM
out_s2_N01000$plot_relbias_dpm
out_s2_N01000$plot_se_dpm
out_s2_N01000$plot_rmse_dpm
out_s2_N01000$plot_power_dpm

# convergence checks
print(out_s2_N01000$flag0_df, n = 100)
print(out_s2_N01000$flag1_df, n = 100)
print(out_s2_N01000$flag2_df, n = 100)

# standard error checks
out_s2_N01000$plot_se_ratio
out_s2_N01000$plot_se_diff

# ML diagnostics
out_s2_N01000$plot_r2
out_s2_N01000$plot_mse


# ---------------------- 00300 -------------------------------

# results CLPM
out_s2_N00300$plot_relbias_clpm
out_s2_N00300$plot_se_clpm
out_s2_N00300$plot_rmse_clpm
out_s2_N00300$plot_power_clpm

# results RI-CLPM
out_s2_N00300$plot_relbias_riclpm
out_s2_N00300$plot_se_riclpm
out_s2_N00300$plot_rmse_riclpm
out_s2_N00300$plot_power_riclpm

# results DPM
out_s2_N00300$plot_relbias_dpm
out_s2_N00300$plot_se_dpm
out_s2_N00300$plot_rmse_dpm
out_s2_N00300$plot_power_dpm

# convergence checks
print(out_s2_N00300$flag0_df, n = 100)
print(out_s2_N00300$flag1_df, n = 100)
print(out_s2_N00300$flag2_df, n = 100)

# standard error checks
out_s2_N00300$plot_se_ratio
out_s2_N00300$plot_se_diff

# ML diagnostics
out_s2_N00300$plot_r2
out_s2_N00300$plot_mse


# ---------------------- 00150 -------------------------------

# results CLPM
out_s2_N00150$plot_relbias_clpm
out_s2_N00150$plot_se_clpm
out_s2_N00150$plot_rmse_clpm
out_s2_N00150$plot_power_clpm

# results RI-CLPM
out_s2_N00150$plot_relbias_riclpm
out_s2_N00150$plot_se_riclpm
out_s2_N00150$plot_rmse_riclpm
out_s2_N00150$plot_power_riclpm

# results DPM
out_s2_N00150$plot_relbias_dpm
out_s2_N00150$plot_se_dpm
out_s2_N00150$plot_rmse_dpm
out_s2_N00150$plot_power_dpm

# convergence checks
print(out_s2_N00150$flag0_df, n = 100)
print(out_s2_N00150$flag1_df, n = 100)
print(out_s2_N00150$flag2_df, n = 100)

# standard error checks
out_s2_N00150$plot_se_ratio
out_s2_N00150$plot_se_diff

# ML diagnostics
out_s2_N00150$plot_r2
out_s2_N00150$plot_mse


# ============================================================
# Scenario 3
# stepwise (rank order stable) non-linear confounding effects
# ============================================================

# ---------------------- 10000 -------------------------------

# results CLPM
out_s3_N10000$plot_relbias_clpm
out_s3_N10000$plot_se_clpm
out_s3_N10000$plot_rmse_clpm
out_s3_N10000$plot_power_clpm

# results RI-CLPM
out_s3_N10000$plot_relbias_riclpm
out_s3_N10000$plot_se_riclpm
out_s3_N10000$plot_rmse_riclpm
out_s3_N10000$plot_power_riclpm

# results DPM
out_s3_N10000$plot_relbias_dpm
out_s3_N10000$plot_se_dpm
out_s3_N10000$plot_rmse_dpm
out_s3_N10000$plot_power_dpm

# convergence checks
print(out_s3_N10000$flag0_df, n = 100)
print(out_s3_N10000$flag1_df, n = 100)
print(out_s3_N10000$flag2_df, n = 100)

# standard error checks
out_s3_N10000$plot_se_ratio
out_s3_N10000$plot_se_diff

# ML diagnostics
out_s3_N10000$plot_r2
out_s3_N10000$plot_mse

# ---------------------- 01000 -------------------------------

# results CLPM
out_s3_N01000$plot_relbias_clpm
out_s3_N01000$plot_se_clpm
out_s3_N01000$plot_rmse_clpm
out_s3_N01000$plot_power_clpm

# results RI-CLPM
out_s3_N01000$plot_relbias_riclpm
out_s3_N01000$plot_se_riclpm
out_s3_N01000$plot_rmse_riclpm
out_s3_N01000$plot_power_riclpm

# results DPM
out_s3_N01000$plot_relbias_dpm
out_s3_N01000$plot_se_dpm
out_s3_N01000$plot_rmse_dpm
out_s3_N01000$plot_power_dpm

# convergence checks
print(out_s3_N01000$flag0_df, n = 100)
print(out_s3_N01000$flag1_df, n = 100)
print(out_s3_N01000$flag2_df, n = 100)

# standard error checks
out_s3_N01000$plot_se_ratio
out_s3_N01000$plot_se_diff

# ML diagnostics
out_s3_N01000$plot_r2
out_s3_N01000$plot_mse


# ---------------------- 00300 -------------------------------

# results CLPM
out_s3_N00300$plot_relbias_clpm
out_s3_N00300$plot_se_clpm
out_s3_N00300$plot_rmse_clpm
out_s3_N00300$plot_power_clpm

# results RI-CLPM
out_s3_N00300$plot_relbias_riclpm
out_s3_N00300$plot_se_riclpm
out_s3_N00300$plot_rmse_riclpm
out_s3_N00300$plot_power_riclpm

# results DPM
out_s3_N00300$plot_relbias_dpm
out_s3_N00300$plot_se_dpm
out_s3_N00300$plot_rmse_dpm
out_s3_N00300$plot_power_dpm

# convergence checks
print(out_s3_N00300$flag0_df, n = 100)
print(out_s3_N00300$flag1_df, n = 100)
print(out_s3_N00300$flag2_df, n = 100)

# standard error checks
out_s3_N00300$plot_se_ratio
out_s3_N00300$plot_se_diff

# ML diagnostics
out_s3_N00300$plot_r2
out_s3_N00300$plot_mse


# ---------------------- 00150 -------------------------------

# results CLPM
out_s3_N00150$plot_relbias_clpm
out_s3_N00150$plot_se_clpm
out_s3_N00150$plot_rmse_clpm
out_s3_N00150$plot_power_clpm

# results RI-CLPM
out_s3_N00150$plot_relbias_riclpm
out_s3_N00150$plot_se_riclpm
out_s3_N00150$plot_rmse_riclpm
out_s3_N00150$plot_power_riclpm

# results DPM
out_s3_N00150$plot_relbias_dpm
out_s3_N00150$plot_se_dpm
out_s3_N00150$plot_rmse_dpm
out_s3_N00150$plot_power_dpm

# convergence checks
print(out_s3_N00150$flag0_df, n = 100)
print(out_s3_N00150$flag1_df, n = 100)
print(out_s3_N00150$flag2_df, n = 100)

# standard error checks
out_s3_N00150$plot_se_ratio
out_s3_N00150$plot_se_diff

# ML diagnostics
out_s3_N00150$plot_r2
out_s3_N00150$plot_mse

# ============================================================
# Scenario 4
# stepwise mixture non-linear confounding effects
# ============================================================

# ---------------------- 10000 -------------------------------

# results CLPM
out_s4_N10000$plot_relbias_clpm
out_s4_N10000$plot_se_clpm
out_s4_N10000$plot_rmse_clpm
out_s4_N10000$plot_power_clpm

# results RI-CLPM
out_s4_N10000$plot_relbias_riclpm
out_s4_N10000$plot_se_riclpm
out_s4_N10000$plot_rmse_riclpm
out_s4_N10000$plot_power_riclpm

# results DPM
out_s4_N10000$plot_relbias_dpm
out_s4_N10000$plot_se_dpm
out_s4_N10000$plot_rmse_dpm
out_s4_N10000$plot_power_dpm

# convergence checks
print(out_s4_N10000$flag0_df, n = 100)
print(out_s4_N10000$flag1_df, n = 100)
print(out_s4_N10000$flag2_df, n = 100)

# standard error checks
out_s4_N10000$plot_se_ratio
out_s4_N10000$plot_se_diff

# ML diagnostics
out_s4_N10000$plot_r2
out_s4_N10000$plot_mse


# ---------------------- 01000 -------------------------------

# results CLPM
out_s4_N01000$plot_relbias_clpm
out_s4_N01000$plot_se_clpm
out_s4_N01000$plot_rmse_clpm
out_s4_N01000$plot_power_clpm

# results RI-CLPM
out_s4_N01000$plot_relbias_riclpm
out_s4_N01000$plot_se_riclpm
out_s4_N01000$plot_rmse_riclpm
out_s4_N01000$plot_power_riclpm

# results DPM
out_s4_N01000$plot_relbias_dpm
out_s4_N01000$plot_se_dpm
out_s4_N01000$plot_rmse_dpm
out_s4_N01000$plot_power_dpm

# convergence checks
print(out_s4_N01000$flag0_df, n = 100)
print(out_s4_N01000$flag1_df, n = 100)
print(out_s4_N01000$flag2_df, n = 100)

# standard error checks
out_s4_N01000$plot_se_ratio
out_s4_N01000$plot_se_diff

# ML diagnostics
out_s4_N01000$plot_r2
out_s4_N01000$plot_mse


# ---------------------- 00300 -------------------------------

# results CLPM
out_s4_N00300$plot_relbias_clpm
out_s4_N00300$plot_se_clpm
out_s4_N00300$plot_rmse_clpm
out_s4_N00300$plot_power_clpm

# results RI-CLPM
out_s4_N00300$plot_relbias_riclpm
out_s4_N00300$plot_se_riclpm
out_s4_N00300$plot_rmse_riclpm
out_s4_N00300$plot_power_riclpm

# results DPM
out_s4_N00300$plot_relbias_dpm
out_s4_N00300$plot_se_dpm
out_s4_N00300$plot_rmse_dpm
out_s4_N00300$plot_power_dpm

# convergence checks
print(out_s4_N00300$flag0_df, n = 100)
print(out_s4_N00300$flag1_df, n = 100)
print(out_s4_N00300$flag2_df, n = 100)

# standard error checks
out_s4_N00300$plot_se_ratio
out_s4_N00300$plot_se_diff

# ML diagnostics
out_s4_N00300$plot_r2
out_s4_N00300$plot_mse


# ---------------------- 00150 -------------------------------

# results CLPM
out_s4_N00150$plot_relbias_clpm
out_s4_N00150$plot_se_clpm
out_s4_N00150$plot_rmse_clpm
out_s4_N00150$plot_power_clpm

# results RI-CLPM
out_s4_N00150$plot_relbias_riclpm
out_s4_N00150$plot_se_riclpm
out_s4_N00150$plot_rmse_riclpm
out_s4_N00150$plot_power_riclpm

# results DPM
out_s4_N00150$plot_relbias_dpm
out_s4_N00150$plot_se_dpm
out_s4_N00150$plot_rmse_dpm
out_s4_N00150$plot_power_dpm

# convergence checks
print(out_s4_N00150$flag0_df, n = 100)
print(out_s4_N00150$flag1_df, n = 100)
print(out_s4_N00150$flag2_df, n = 100)

# standard error checks
out_s4_N00150$plot_se_ratio
out_s4_N00150$plot_se_diff

# ML diagnostics
out_s4_N00150$plot_r2
out_s4_N00150$plot_mse


# ============================================================
# Scenario 5
# stepwise (rank order stable) non-linear confounding effects
# omitting c1 and c2
# ============================================================

# ---------------------- 10000 -------------------------------

# results CLPM
out_s5_N10000$plot_relbias_clpm
out_s5_N10000$plot_se_clpm
out_s5_N10000$plot_rmse_clpm
out_s5_N10000$plot_power_clpm

# results RI-CLPM
out_s5_N10000$plot_relbias_riclpm
out_s5_N10000$plot_se_riclpm
out_s5_N10000$plot_rmse_riclpm
out_s5_N10000$plot_power_riclpm

# results DPM
out_s5_N10000$plot_relbias_dpm
out_s5_N10000$plot_se_dpm
out_s5_N10000$plot_rmse_dpm
out_s5_N10000$plot_power_dpm

# convergence checks
print(out_s5_N10000$flag0_df, n = 100)
print(out_s5_N10000$flag1_df, n = 100)
print(out_s5_N10000$flag2_df, n = 100)

# standard error checks
out_s5_N10000$plot_se_ratio
out_s5_N10000$plot_se_diff

# ML diagnostics
out_s5_N10000$plot_r2
out_s5_N10000$plot_mse


# ---------------------- 01000 -------------------------------

# results CLPM
out_s5_N01000$plot_relbias_clpm
out_s5_N01000$plot_se_clpm
out_s5_N01000$plot_rmse_clpm
out_s5_N01000$plot_power_clpm
 
# results RI-CLPM
out_s5_N01000$plot_relbias_riclpm
out_s5_N01000$plot_se_riclpm
out_s5_N01000$plot_rmse_riclpm
out_s5_N01000$plot_power_riclpm

# results DPM
out_s5_N01000$plot_relbias_dpm
out_s5_N01000$plot_se_dpm
out_s5_N01000$plot_rmse_dpm
out_s5_N01000$plot_power_dpm

# convergence checks
print(out_s5_N01000$flag0_df, n = 100)
print(out_s5_N01000$flag1_df, n = 100)
print(out_s5_N01000$flag2_df, n = 100)

# standard error checks
out_s5_N01000$plot_se_ratio
out_s5_N01000$plot_se_diff

# ML diagnostics
out_s5_N01000$plot_r2
out_s5_N01000$plot_mse


# ---------------------- 00300 -------------------------------

# results CLPM
out_s5_N00300$plot_relbias_clpm
out_s5_N00300$plot_se_clpm
out_s5_N00300$plot_rmse_clpm
out_s5_N00300$plot_power_clpm

# results RI-CLPM
out_s5_N00300$plot_relbias_riclpm
out_s5_N00300$plot_se_riclpm
out_s5_N00300$plot_rmse_riclpm
out_s5_N00300$plot_power_riclpm

# results DPM
out_s5_N00300$plot_relbias_dpm
out_s5_N00300$plot_se_dpm
out_s5_N00300$plot_rmse_dpm
out_s5_N00300$plot_power_dpm

# convergence checks
print(out_s5_N00300$flag0_df, n = 100)
print(out_s5_N00300$flag1_df, n = 100)
print(out_s5_N00300$flag2_df, n = 100)

# standard error checks
out_s5_N00300$plot_se_ratio
out_s5_N00300$plot_se_diff

# ML diagnostics
out_s5_N00300$plot_r2
out_s5_N00300$plot_mse


# ---------------------- 00150 -------------------------------

# results CLPM
out_s5_N00150$plot_relbias_clpm
out_s5_N00150$plot_se_clpm
out_s5_N00150$plot_rmse_clpm
out_s5_N00150$plot_power_clpm

# results RI-CLPM
out_s5_N00150$plot_relbias_riclpm
out_s5_N00150$plot_se_riclpm
out_s5_N00150$plot_rmse_riclpm
out_s5_N00150$plot_power_riclpm

# results DPM
out_s5_N00150$plot_relbias_dpm
out_s5_N00150$plot_se_dpm
out_s5_N00150$plot_rmse_dpm
out_s5_N00150$plot_power_dpm

# convergence checks
print(out_s5_N00150$flag0_df, n = 100)
print(out_s5_N00150$flag1_df, n = 100)
print(out_s5_N00150$flag2_df, n = 100)

# standard error checks
out_s5_N00150$plot_se_ratio
out_s5_N00150$plot_se_diff

# ML diagnostics
out_s5_N00150$plot_r2
out_s5_N00150$plot_mse


# ============================================================
# Scenario 6
# stepwise mixture non-linear confounding effects
# omitting c1 and c2
# ============================================================


# ---------------------- 10000 -------------------------------

# results CLPM
out_s6_N10000$plot_relbias_clpm
out_s6_N10000$plot_se_clpm
out_s6_N10000$plot_rmse_clpm
out_s6_N10000$plot_power_clpm

# results RI-CLPM
out_s6_N10000$plot_relbias_riclpm
out_s6_N10000$plot_se_riclpm
out_s6_N10000$plot_rmse_riclpm
out_s6_N10000$plot_power_riclpm

# results DPM
out_s6_N10000$plot_relbias_dpm
out_s6_N10000$plot_se_dpm
out_s6_N10000$plot_rmse_dpm
out_s6_N10000$plot_power_dpm

# convergence checks
print(out_s6_N10000$flag0_df, n = 100)
print(out_s6_N10000$flag1_df, n = 100)
print(out_s6_N10000$flag2_df, n = 100)

# standard error checks
out_s6_N10000$plot_se_ratio
out_s6_N10000$plot_se_diff

# ML diagnostics
out_s6_N10000$plot_r2
out_s6_N10000$plot_mse


# ---------------------- 01000 -------------------------------

# results CLPM
out_s6_N01000$plot_relbias_clpm
out_s6_N01000$plot_se_clpm
out_s6_N01000$plot_rmse_clpm
out_s6_N01000$plot_power_clpm

# results RI-CLPM
out_s6_N01000$plot_relbias_riclpm
out_s6_N01000$plot_se_riclpm
out_s6_N01000$plot_rmse_riclpm
out_s6_N01000$plot_power_riclpm

# results DPM
out_s6_N01000$plot_relbias_dpm
out_s6_N01000$plot_se_dpm
out_s6_N01000$plot_rmse_dpm
out_s6_N01000$plot_power_dpm

# convergence checks
print(out_s6_N01000$flag0_df, n = 100)
print(out_s6_N01000$flag1_df, n = 100)
print(out_s6_N01000$flag2_df, n = 100)

# standard error checks
out_s6_N01000$plot_se_ratio
out_s6_N01000$plot_se_diff

# ML diagnostics
out_s6_N01000$plot_r2
out_s6_N01000$plot_mse


# ---------------------- 00300 -------------------------------

# results CLPM
out_s6_N00300$plot_relbias_clpm
out_s6_N00300$plot_se_clpm
out_s6_N00300$plot_rmse_clpm
out_s6_N00300$plot_power_clpm

# results RI-CLPM
out_s6_N00300$plot_relbias_riclpm
out_s6_N00300$plot_se_riclpm
out_s6_N00300$plot_rmse_riclpm
out_s6_N00300$plot_power_riclpm

# results DPM
out_s6_N00300$plot_relbias_dpm
out_s6_N00300$plot_se_dpm
out_s6_N00300$plot_rmse_dpm
out_s6_N00300$plot_power_dpm

# convergence checks
print(out_s6_N00300$flag0_df, n = 100)
print(out_s6_N00300$flag1_df, n = 100)
print(out_s6_N00300$flag2_df, n = 100)

# standard error checks
out_s6_N00300$plot_se_ratio
out_s6_N00300$plot_se_diff

# ML diagnostics
out_s6_N00300$plot_r2
out_s6_N00300$plot_mse


# ---------------------- 00150 -------------------------------

# results CLPM
out_s6_N00150$plot_relbias_clpm
out_s6_N00150$plot_se_clpm
out_s6_N00150$plot_rmse_clpm
out_s6_N00150$plot_power_clpm

# results RI-CLPM
out_s6_N00150$plot_relbias_riclpm
out_s6_N00150$plot_se_riclpm
out_s6_N00150$plot_rmse_riclpm
out_s6_N00150$plot_power_riclpm

# results DPM
out_s6_N00150$plot_relbias_dpm
out_s6_N00150$plot_se_dpm
out_s6_N00150$plot_rmse_dpm
out_s6_N00150$plot_power_dpm

# convergence checks
print(out_s6_N00150$flag0_df, n = 100)
print(out_s6_N00150$flag1_df, n = 100)
print(out_s6_N00150$flag2_df, n = 100)

# standard error checks
out_s6_N00150$plot_se_ratio
out_s6_N00150$plot_se_diff

# ML diagnostics
out_s6_N00150$plot_r2
out_s6_N00150$plot_mse

