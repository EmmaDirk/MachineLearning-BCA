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
source(here("01_scripts", "02_development", "06_thesis_v3", "12_plotting.R"))

# ============================================================
# load data
# ============================================================

# N = 10000
dat_N10000_s1 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s01_N10000.rds"))
dat_N10000_s2 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s02_N10000.rds"))
dat_N10000_s3 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s03_N10000.rds"))
dat_N10000_s4 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s04_N10000.rds"))
dat_N10000_s5 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s05_N10000.rds"))
dat_N10000_s6 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s06_N10000.rds"))

# N = 1000
dat_N01000_s1 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s01_N01000.rds"))
dat_N01000_s2 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s02_N01000.rds"))
dat_N01000_s3 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s03_N01000.rds"))
dat_N01000_s4 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s04_N01000.rds"))
dat_N01000_s5 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s05_N01000.rds"))
dat_N01000_s6 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s06_N01000.rds"))

# N = 300
dat_N00300_s1 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s01_N00300.rds"))
dat_N00300_s2 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s02_N00300.rds"))
dat_N00300_s3 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s03_N00300.rds"))
dat_N00300_s4 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s04_N00300.rds"))
dat_N00300_s5 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s05_N00300.rds"))
dat_N00300_s6 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s06_N00300.rds"))

# N = 150
dat_N00150_s1 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s01_N00150.rds"))
dat_N00150_s2 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s02_N00150.rds"))
dat_N00150_s3 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s03_N00150.rds"))
dat_N00150_s4 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s04_N00150.rds"))
dat_N00150_s5 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s05_N00150.rds"))
dat_N00150_s6 <- readRDS(here("02_data", "02_server_runs", "02_run", "test_s06_N00150.rds"))

# ============================================================
# results
# ============================================================
out_s1_N10000 <- plot_overview_suite(dat_N10000_s1)
out_s1_N01000 <- plot_overview_suite(dat_N01000_s1)
out_s1_N00300 <- plot_overview_suite(dat_N00300_s1)
out_s1_N00150 <- plot_overview_suite(dat_N00150_s1)

out_s2_N10000 <- plot_overview_suite(dat_N10000_s2)
out_s2_N01000 <- plot_overview_suite(dat_N01000_s2)
out_s2_N00300 <- plot_overview_suite(dat_N00300_s2)
out_s2_N00150 <- plot_overview_suite(dat_N00150_s2)

out_s3_N10000 <- plot_overview_suite(dat_N10000_s3)
out_s3_N01000 <- plot_overview_suite(dat_N01000_s3)
out_s3_N00300 <- plot_overview_suite(dat_N00300_s3)
out_s3_N00150 <- plot_overview_suite(dat_N00150_s3)

out_s4_N10000 <- plot_overview_suite(dat_N10000_s4)
out_s4_N01000 <- plot_overview_suite(dat_N01000_s4)
out_s4_N00300 <- plot_overview_suite(dat_N00300_s4)
out_s4_N00150 <- plot_overview_suite(dat_N00150_s4)

out_s5_N10000 <- plot_overview_suite(dat_N10000_s5)
out_s5_N01000 <- plot_overview_suite(dat_N01000_s5)
out_s5_N00300 <- plot_overview_suite(dat_N00300_s5)
out_s5_N00150 <- plot_overview_suite(dat_N00150_s5)

out_s6_N10000 <- plot_overview_suite(dat_N10000_s6)
out_s6_N01000 <- plot_overview_suite(dat_N01000_s6)
out_s6_N00300 <- plot_overview_suite(dat_N00300_s6)
out_s6_N00150 <- plot_overview_suite(dat_N00150_s6)

# 
