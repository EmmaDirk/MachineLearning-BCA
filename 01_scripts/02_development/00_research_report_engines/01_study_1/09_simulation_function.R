# this script runs the wrapper function that executes one replication of the simulation study
# adding a progress bar and parallelization
# ------------------------------------------------------------------------------------------------

run_simulation_study <- function(
  reps,                                                                    # number of replications
  N,                                                                       # sample size
  T,                                                                       # number of waves
  k,                                                                       # number of confounders
  scenarios,                                                               # e.g., c("constant","stepwise")
  B_scenarios,                                                             # named list of beta trajectories
  A,                                                                       # 2×2 AR + cross-lag matrix
  Psi,                                                                     # k×k confounder covariance
  rho_extra,                                                               # extra covariance added to observations
  models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca")
  cores = NULL,                                                            # default is detectCores()/2
  base_seed = 1234                                                         # master seed for reproducible reps
) {

  # choose number of cores
  if (is.null(cores)) {
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # run sequentially if cores is 1
  if (cores == 1L) {

    results_list <- lapply(
      X = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study(
          rep_id        = rep_id,
          N             = N,
          T             = T,
          k             = k,
          scenarios     = scenarios,
          B_scenarios   = B_scenarios,
          A             = A,
          Psi           = Psi,
          rho_extra     = rho_extra,
          models_to_run = models_to_run,
          base_seed     = base_seed
        )
      }
    )

    return(dplyr::bind_rows(results_list))
  }

  # make the cluster
  cl <- parallel::makeCluster(cores)

  # set RNG streams on the cluster
  parallel::clusterSetRNGStream(cl, iseed = base_seed)

  # load packages on each worker
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(mvtnorm)
    NULL
  })

  # export functions and objects needed on workers
  parallel::clusterExport(
    cl,
    c(
      "simulate_panel_data",
      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",
      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "extract_lagged_parameters",
      "extract_rho_vec",
      "residualise_panel_linearC",
      "run_one_rep_study",
      "N","T","k","scenarios","B_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study(
        rep_id        = rep_id,
        N             = N,
        T             = T,
        k             = k,
        scenarios     = scenarios,
        B_scenarios   = B_scenarios,
        A             = A,
        Psi           = Psi,
        rho_extra     = rho_extra,
        models_to_run = models_to_run,
        base_seed     = base_seed
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # bind and return
  dplyr::bind_rows(results_list)
}
