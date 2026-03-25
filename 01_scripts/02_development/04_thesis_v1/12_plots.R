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

# turn into one long dataframe
all_model_dfs_300_long <- bind_rows(all_model_300_dfs, .id = "model")

# plot the results
plot_engine_results(all_model_dfs_long)
