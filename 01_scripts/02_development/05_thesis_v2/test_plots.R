# ============================================================
# libraries + setup
# ============================================================
library(here)
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))

# source plotting function
source(here("01_scripts", "02_development", "04_thesis_v1", "11_plotting.R"))

# ============================================================
# load data
# ============================================================
all_model_300_dfs <- readRDS(
  here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_300.rds"
  )
)

all_model_5000_dfs <- readRDS(
  here::here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_5000_model_set.rds"
  )
)

all_model_300_df <- dplyr::bind_rows(all_model_300_dfs)
all_model_5000_df <- dplyr::bind_rows(all_model_5000_dfs)

# ============================================================
# results
# ============================================================
out_300 <- plot_engine_results(all_model_300_df)
out_5000 <- plot_engine_results(all_model_5000_df)

# 1) plots for n = 300
out_300$plot_relbias
out_300$plot_rmse
out_300$plot_se

out_300$plot_relbias_clpm
out_300$plot_rmse_clpm
out_300$plot_se_clpm

out_300$plot_relbias_riclpm
out_300$plot_rmse_riclpm
out_300$plot_se_riclpm

out_300$plot_relbias_dpm
out_300$plot_rmse_dpm
out_300$plot_se_dpm

# 2) plots for n = 5000
out_5000$plot_relbias
out_5000$plot_rmse
out_5000$plot_se

out_5000$plot_relbias_clpm
out_5000$plot_rmse_clpm
out_5000$plot_se_clpm

out_5000$plot_relbias_riclpm
out_5000$plot_rmse_riclpm
out_5000$plot_se_riclpm

out_5000$plot_relbias_dpm
out_5000$plot_rmse_dpm
out_5000$plot_se_dpm

# ============================================================
# diagnostics
# ============================================================

# 1) plots for n = 300
out_300$plot_se_ratio
out_300$plot_se_diff
out_300$plot_flag0
out_300$plot_flag1
out_300$plot_flag2
out_300$plot_bootstrap_prop_success
out_300$plot_improper_reason_top
out_300$plot_bootstrap_issue_top

# 2) plots for n = 5000
out_5000$plot_se_ratio
out_5000$plot_se_diff
out_5000$plot_flag0
out_5000$plot_flag1
out_5000$plot_flag2
out_5000$plot_bootstrap_prop_success
out_5000$plot_improper_reason_top
out_5000$plot_bootstrap_issue_top

# ============================================================
# tables / summaries
# ============================================================

# 1) summaries for n = 300
out_300$flag_summary_method
out_300$flag_summary_method_T
out_300$improper_reason_df
out_300$improper_reason_top_df
out_300$bootstrap_issue_df
out_300$bootstrap_issue_top_df
out_300$se_check
out_300$method_key

# 2) summaries for n = 5000
out_5000$flag_summary_method
out_5000$flag_summary_method_T
out_5000$improper_reason_df
out_5000$improper_reason_top_df
out_5000$bootstrap_issue_df
out_5000$bootstrap_issue_top_df
out_5000$se_check
out_5000$method_key
