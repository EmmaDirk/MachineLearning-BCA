# script to test if the lavaan models that are built are working. 
# there are 6 models:
# - CLPM without confounder adjustment
# - CLPM with direct confounder adjustment
# - RI-CLPM with indirect confounder adjustment via random intercepts
# - DPM
# - RI-CLPM with freed latent factor loadings
# - DPM with freed latent factor loadings
# ------------------------------------------------------------------------------------------------------------

library(here)
library(lavaan)

# load the model strings
source(here("01_scripts", "03_exploration", "02_sim_engine", "04_models", "v1_lavaan.R"))

# build a CLPM model string without confounder adjustment
clpm <- build_clpm(T = 5, k=0, confounder_order = 0)
cat(clpm)

# build a clpm that controls only for main effects 
clpm <- build_clpm(T = 5, k=3, confounder_order = 1)
cat(clpm)

# build a clpm that controls for main effects and 2-way interactions
clpm <- build_clpm(T = 5, k=3, confounder_order = 2)
cat(clpm)

# build a clpm that controls for main effects and 2-way + 3-way interactions
clpm <- build_clpm(T = 5, k=3, confounder_order = 3)
cat(clpm)

# build a CLPM with main effects, but exclude c1
clpm <- build_clpm(T = 5, k=3, confounder_order = 1, exclude = "c1")
cat(clpm)

# build a RI-CLPM
ri_clpm <- build_riclpm(T = 5, free_loadings = FALSE)
cat(ri_clpm)

# build a RI-CLPM with freed loadings
ri_clpm <- build_riclpm(T = 5, free_loadings = TRUE)
cat(ri_clpm)

# build a DPM
dpm <- build_dpm(T = 5, free_loadings = FALSE)
cat(dpm)

# build a DPM with freed loadings
dpm <- build_dpm(T = 5, free_loadings = TRUE)
cat(dpm)
