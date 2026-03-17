# This script sources the engine scripts in the order they depend on each other.
# The purpose is simply convenience.
# If you prefer, you can also source the files one by one.
# ------------------------------------------------------------------------------------------

source("00_packages.R")
source("01_delta_sampler.R")
source("02_delta_trajectory.R")
source("03_data_generation.R")
source("04_model_strings.R")
source("05_residualisers.R")
source("06_fit_helpers.R")
source("07_fitters.R")
source("08_simulation_helpers.R")
source("09_simulation_runner.R")
