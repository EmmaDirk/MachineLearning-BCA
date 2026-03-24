# library
library(here)
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))

# source the plotting function
source(here("01_scripts", "02_development", "04_thesis_v1", "11_plotting.R"))

# first load the data we want to use:
dat_00 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "00_none_CLPMno_constant_1c_linear.rds"))
dat_01 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "01_BCALM_CLPM_constant_1c_linear.rds"))
dat_02 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "02_none_CLPM_constant_1c_linear.rds"))
dat_03 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "03_BCAXGB_CLPM_constant_1c_linear.rds"))
dat_04 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "04_none_RICLPM_constant_1c_linear.rds"))
dat_05 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "05_none_DPM_constant_1c_linear.rds"))

# process data
dat_00 <- dat_00$results
dat_01 <- dat_01$results
dat_02 <- dat_02$results
dat_03 <- dat_03$results
dat_04 <- dat_04$results
dat_05 <- dat_05$results

# merge
dat_all <- rbind(dat_00, dat_01, dat_02, dat_03, dat_04, dat_05)

# plot the data
plot_engine_results(dat_all)

# repeat the story for N=300
dat_003 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "00_none_CLPMno_constant_1c_linear_300.rds"))
dat_013 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "01_BCALM_CLPM_constant_1c_linear_300.rds"))
dat_023 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "02_none_CLPM_constant_1c_linear_300.rds"))
dat_033 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "03_BCAXGB_CLPM_constant_1c_linear_300.rds"))
dat_043 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "04_none_RICLPM_constant_1c_linear_300.rds"))
dat_053 <- readRDS(here("03_output", "02_thesis", "01_tests", "01_data", "05_none_DPM_constant_1c_linear_300.rds"))

# process data
dat_003 <- dat_003$results
dat_013 <- dat_013$results
dat_023 <- dat_023$results
dat_033 <- dat_033$results
dat_043 <- dat_043$results
dat_053 <- dat_053$results

# merge
dat_all3 <- rbind(dat_003, dat_013, dat_023, dat_033, dat_043, dat_053)

# plot the data
plot_engine_results(dat_all3)