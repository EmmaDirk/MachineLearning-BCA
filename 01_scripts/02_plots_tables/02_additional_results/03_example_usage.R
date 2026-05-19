# =================================================================================================
#
# This script gives a compact example of how to inspect additional simulation results.
#
# It loads the combined simulation-results data frame, sources the overview functions, creates one
# overview object for a selected scenario and sample size, and displays the tables and plots returned
# by that object.
#
# The example is intended as a quick entry point for browsing the available additional results.
# =================================================================================================

# ---- setup ---------------------------------------------------------------------------------------

library(here)
library(tidyverse)
library(scales)

# Define the project root.

study_root <- here::here()

# Define reusable directories.

script_dir <- file.path(
  study_root,
  "01_scripts",
  "02_plots_tables",
  "02_additional_results"
)

data_dir <- file.path(
  study_root,
  "02_data",
  "02_additional_results"
)

# Source the overview functions.

source(file.path(script_dir, "01_overview_models.R"))
source(file.path(script_dir, "02_plot_overview_suite.R"))


# ---- load data -----------------------------------------------------------------------------------

# Load the combined simulation-results data frame.

dat_all <- readRDS(file.path(data_dir, "s1234_N300_1000_2000.rds"))


# ---- inspect model specifications ---------------------------------------------------------------

# Show the model specifications present in the combined results file.

models_all <- overview_models(
  results_df = dat_all,
  include_counts = TRUE,
  print_result = TRUE
)


# ---- create overview object ---------------------------------------------------------------------

# Create the overview object for one scenario and one sample size.

out_s1_N1000 <- plot_overview_suite(
  results_df = dat_all,
  scenario_id = 1,
  N = 1000,
  occasions = 2:5,
  drop_flagged = TRUE
)


# ---- returned tables ----------------------------------------------------------------------------

# Model specifications included in the overview object.

out_s1_N1000$model_overview

# Relative bias for each SEM parameter, method, model family, and occasion.

out_s1_N1000$relbias_df

# Root mean squared error for each SEM parameter, method, model family, and occasion.

out_s1_N1000$rmse_df

# Mean reported standard error for each SEM parameter, method, model family, and occasion.

out_s1_N1000$se_df

# Detection probability for each SEM parameter, method, model family, and occasion.

out_s1_N1000$detect_df

# Standard-error calibration table comparing empirical Monte Carlo standard deviation with mean
# reported standard error.

out_s1_N1000$se_check

# Prediction MSE for residualizer diagnostics.

out_s1_N1000$mse_df

# Prediction R-squared for residualizer diagnostics.

out_s1_N1000$r2_df

# Theoretical reference values for the residualizer diagnostics.

out_s1_N1000$theory_df

# True confounder-effect trajectories extracted from true_delta_t_vector.

out_s1_N1000$true_delta_df

# Proportion of flag0 by method.

out_s1_N1000$flag0_df

# Proportion of flag1 by method.

out_s1_N1000$flag1_df

# Proportion of flag2 by method.

out_s1_N1000$flag2_df

# Breakdown of improper-solution reasons by method.

out_s1_N1000$improper_reason_df


# ---- CLPM plots ----------------------------------------------------------------------------------

# Relative bias for CLPM-family methods.

out_s1_N1000$plot_relbias_clpm

# Mean reported standard errors for CLPM-family methods.

out_s1_N1000$plot_se_clpm

# Root mean squared error for CLPM-family methods.

out_s1_N1000$plot_rmse_clpm

# Detection probability for CLPM-family methods.

out_s1_N1000$plot_power_clpm


# ---- RI-CLPM plots -------------------------------------------------------------------------------

# Relative bias for RI-CLPM-family methods.

out_s1_N1000$plot_relbias_riclpm

# Mean reported standard errors for RI-CLPM-family methods.

out_s1_N1000$plot_se_riclpm

# Root mean squared error for RI-CLPM-family methods.

out_s1_N1000$plot_rmse_riclpm

# Detection probability for RI-CLPM-family methods.

out_s1_N1000$plot_power_riclpm


# ---- DPM plots -----------------------------------------------------------------------------------

# Relative bias for DPM-family methods.

out_s1_N1000$plot_relbias_dpm

# Mean reported standard errors for DPM-family methods.

out_s1_N1000$plot_se_dpm

# Root mean squared error for DPM-family methods.

out_s1_N1000$plot_rmse_dpm

# Detection probability for DPM-family methods.

out_s1_N1000$plot_power_dpm


# ---- standard-error calibration plots -----------------------------------------------------------

# Ratio of mean reported standard error to empirical Monte Carlo standard deviation for CLPM-family
# methods.

out_s1_N1000$plot_se_ratio_clpm

# Ratio of mean reported standard error to empirical Monte Carlo standard deviation for RI-CLPM-family
# methods.

out_s1_N1000$plot_se_ratio_riclpm

# Ratio of mean reported standard error to empirical Monte Carlo standard deviation for DPM-family
# methods.

out_s1_N1000$plot_se_ratio_dpm

# Difference between mean reported standard error and empirical Monte Carlo standard deviation.

out_s1_N1000$plot_se_diff


# ---- residualizer diagnostic plots --------------------------------------------------------------

# Prediction R-squared for the residualizer part of the analysis.

out_s1_N1000$plot_r2

# Prediction MSE for the residualizer part of the analysis.

out_s1_N1000$plot_mse


# ---- true-delta plot -----------------------------------------------------------------------------

# True confounder-effect trajectories over time.

out_s1_N1000$plot_true_delta


# ---- flag plots ----------------------------------------------------------------------------------

# Proportion of flag0 by method.

out_s1_N1000$plot_flag0

# Proportion of flag1 by method.

out_s1_N1000$plot_flag1

# Proportion of flag2 by method.

out_s1_N1000$plot_flag2

# Improper-solution reasons by method.

out_s1_N1000$plot_improper_reasons


# ---- returned object names -----------------------------------------------------------------------

# Names of all top-level objects returned by plot_overview_suite().

names(out_s1_N1000)