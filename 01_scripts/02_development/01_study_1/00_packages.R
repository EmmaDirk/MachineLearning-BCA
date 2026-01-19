# this script loads the required packages for this simulation study
# we do not install packages since they are managed by renv
# ---------------------------------------------------------------------

# required packages
pkgs <- c("mvtnorm", "lavaan", "tidyverse", "here", 
          "parallel", "pbapply", "ggh4x", 
          "colorspace", "patchwork", "viridis")

# load packages
lapply(pkgs, library, character.only = TRUE)
