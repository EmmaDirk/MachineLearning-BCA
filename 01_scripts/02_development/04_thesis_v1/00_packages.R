# This script contains very small helpers for package management.
# The goal is not to hide anything. The goal is simply to keep the package checks
# in one obvious place, so the other scripts can stay focused on their real job.
#
# We do two things here:
# 1) check whether the needed packages are installed
# 2) optionally attach the packages to the search path
#
# The engine itself mostly uses requireNamespace() inside the relevant functions.
# That means you do not strictly need to call these helpers.
# Still, for a clean session, it is often useful to call them once at the start.
# ------------------------------------------------------------------------------------------

# helper that checks whether a package is installed
is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# check the packages that the engine may need
check_sim_engine_packages <- function(include_xgboost = FALSE) {

  # core packages used by the simulation and SEM parts
  needed <- c("mvtnorm", "lavaan", "pbapply")

  # xgboost is optional, because only the xgb residualiser needs it
  if (isTRUE(include_xgboost)) {
    needed <- c(needed, "xgboost")
  }

  # find which packages are missing
  missing <- needed[!vapply(needed, is_installed, logical(1))]

  # return a small summary instead of stopping immediately
  list(
    needed = needed,
    missing = missing,
    all_available = length(missing) == 0
  )
}

# optionally attach the packages for interactive work
load_sim_engine_packages <- function(include_xgboost = FALSE) {

  chk <- check_sim_engine_packages(include_xgboost = include_xgboost)

  # if something is missing, stop with a clear message
  if (!chk$all_available) {
    stop(
      "The following package(s) are missing: ",
      paste(chk$missing, collapse = ", "),
      ". Please install them first."
    )
  }

  # attach the core packages
  library(mvtnorm)
  library(lavaan)
  library(pbapply)

  # attach xgboost only if requested
  if (isTRUE(include_xgboost)) {
    library(xgboost)
  }

  invisible(chk)
}
