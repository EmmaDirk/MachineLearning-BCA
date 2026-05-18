# ---------------------------------------------------------------------
# this script loads the required packages for this simulation study
# ---------------------------------------------------------------------

# required packages
pkgs <- c("mvtnorm", "lavaan", "tidyverse", "here", 
          "parallel", "pbapply", "ggh4x", 
          "patchwork", "xgboost", "glmnet")

# load packages
lapply(pkgs, library, character.only = TRUE)
