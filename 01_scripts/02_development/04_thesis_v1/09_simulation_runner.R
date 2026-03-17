# This script runs the full simulation study by repeatedly calling the direct
# model-fitting functions through a small specification list.
#
# The key difference from the older setup is this:
# - the model functions themselves now do the real work
# - the simulation runner only repeats those calls over replications
# - this keeps the simulation code relatively light
#
# Inputs:
# - reps            : number of replications
# - N               : sample size
# - T               : number of waves
# - Delta_scenarios : named list of delta trajectories
# - Phi             : 2 x 2 lag matrix
# - Omega11         : covariance matrix of the base confounders
# - Sigma           : target covariance matrix of (X_t, Y_t)
# - model_specs     : named list describing which model functions to call
#
# Output:
# - one tidy data frame with one row per estimated lagged parameter
# ------------------------------------------------------------------------------------------

run_simulation_study <- function(
  reps,
  N,
  T,
  Delta_scenarios,
  Phi,
  Omega11,
  Sigma,
  model_specs,
  cores = NULL,
  base_seed = 1234
) {

  # ------------------------- input checks -------------------------

  if (!is_whole_number_scalar(reps) || reps < 1) {
    stop("'reps' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(N) || N < 1) {
    stop("'N' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(T) || T < 2) {
    stop("'T' must be a single integer >= 2.")
  }

  if (!is.list(Delta_scenarios) || length(Delta_scenarios) == 0) {
    stop("'Delta_scenarios' must be a non-empty named list.")
  }

  if (is.null(cores)) {
    cores <- max(1L, floor(parallel::detectCores() / 2))
  }

  if (!is_whole_number_scalar(cores) || cores < 1) {
    stop("'cores' must be a single integer >= 1.")
  }

  # ------------------------- sequential case -------------------------
  # when cores = 1, the function stays very simple.

  if (cores == 1L) {

    results_list <- pbapply::pblapply(
      X = seq_len(reps),
      FUN = function(rep_id) {
        run_one_replication(
          rep_id = rep_id,
          N = N,
          T = T,
          Delta_scenarios = Delta_scenarios,
          Phi = Phi,
          Omega11 = Omega11,
          Sigma = Sigma,
          model_specs = model_specs,
          base_seed = base_seed
        )
      }
    )

    return(do.call(rbind, results_list))
  }

  # ------------------------- parallel case -------------------------
  # we use a PSOCK cluster.
  # each worker receives the functions it needs and then calls run_one_replication().

  cl <- parallel::makeCluster(cores)

  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterSetRNGStream(cl, iseed = base_seed)

  parallel::clusterEvalQ(cl, {
    library(mvtnorm)
    library(lavaan)
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      simulation_engine_exports(),
      "N",
      "T",
      "Delta_scenarios",
      "Phi",
      "Omega11",
      "Sigma",
      "model_specs",
      "base_seed"
    ),
    envir = environment()
  )

  results_list <- pbapply::pblapply(
    X = seq_len(reps),
    cl = cl,
    FUN = function(rep_id) {
      run_one_replication(
        rep_id = rep_id,
        N = N,
        T = T,
        Delta_scenarios = Delta_scenarios,
        Phi = Phi,
        Omega11 = Omega11,
        Sigma = Sigma,
        model_specs = model_specs,
        base_seed = base_seed
      )
    }
  )

  do.call(rbind, results_list)
}
