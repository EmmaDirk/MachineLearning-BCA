# library
library(here)
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))

# source the plotting function
source(here("01_scripts", "02_development", "04_thesis_v1", "11_plotting.R"))

# get the data
all_model_300_dfs <- readRDS(
  here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_300.rds"
  )
)

all_model_5000_dfs <- readRDS(
  here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_5000.rds"
  )
)

# transform to the long format
all_model_300_df <- dplyr::bind_rows(all_model_300_dfs)
all_model_5000_df <- dplyr::bind_rows(all_model_5000_dfs)

# plot the results
out_300 <- plot_engine_results(all_model_300_df)
out_5000 <- plot_engine_results(all_model_5000_df)

# 1) show the plots for sample size 300
# a) complete plots
out_300$plot_relbias
out_300$plot_rmse
out_300$plot_se

# b) clpm plots
out_300$plot_relbias_clpm
out_300$plot_rmse_clpm
out_300$plot_se_clpm

# c) riclpm plots
out_300$plot_relbias_riclpm
out_300$plot_rmse_riclpm
out_300$plot_se_riclpm

# d) dpm plots
out_300$plot_relbias_dpm
out_300$plot_rmse_dpm
out_300$plot_se_dpm

# 2) show the plots for sample size 5000
# a) complete plots
out_5000$plot_relbias
out_5000$plot_rmse
out_5000$plot_se

# b) clpm plots
out_5000$plot_relbias_clpm
out_5000$plot_rmse_clpm
out_5000$plot_se_clpm

# c) riclpm plots
out_5000$plot_relbias_riclpm
out_5000$plot_rmse_riclpm
out_5000$plot_se_riclpm

# d) dpm plots
out_5000$plot_relbias_dpm
out_5000$plot_rmse_dpm
out_5000$plot_se_dpm
