# This function runs the full simulation study by repeatedly calling run_one_rep_study() across replications,
# either sequentially (single core) or in parallel (multiple cores).
#
# 1) It determines how many CPU cores to use:
# - if cores is NULL, it defaults to detectCores()/2 (at least 1)
# - if cores == 1, it runs everything sequentially via lapply()
#
# 2) If running sequentially (cores == 1):
# - for rep_id = 1..reps, it calls run_one_rep_study() with the provided settings
# - it row-binds the per-replication outputs into one combined results data frame
#
# 3) If running in parallel (cores > 1):
# - it creates a PSOCK cluster with the requested number of workers
# - it initializes reproducible random number streams on the cluster (clusterSetRNGStream with base_seed)
# - it loads required packages on each worker (lavaan, mvtnorm)
# - it exports all functions and objects needed by the workers (including D_scenarios and model settings)
#
# 4) It distributes replications across workers:
# - it uses pbapply::pblapply() to run rep_id = 1..reps in parallel with a progress bar
# - each worker calls run_one_rep_study(), which handles scenario checks, data simulation, model fitting, and extraction
#
# 5) It finalizes and returns results:
# - it stops the cluster to free resources
# - it row-binds all replication outputs into a single results data frame and returns it
# ------------------------------------------------------------------------------------------------------------

run_simulation_study <- function(
  reps,                                                                    # number of replications
  N,                                                                       # sample size
  T,                                                                       # number of waves
  k,                                                                       # number of confounders
  scenarios,                                                               # e.g., c("constant","stepwise")
  D_scenarios,                                                             # named list of delta trajectories
  A,                                                                       # 2×2 autoregressive (beta) + cross-lagged (gamma) matrix
  Psi,                                                                     # k×k confounder covariance
  rho_extra,                                                               # extra covariance added to observations
  models_to_run,                                                           # c("clpm","riclpm","riclpm_free_RI_loadings","dpm","dpm_free_loadings","adj","lbca")
  cores = NULL,                                                            # default is detectCores()/2
  base_seed = 1234,                                                        # master seed
  ci_level = 0.95                                                          # CI level for extracted parameters
){

  # choose number of cores
  if (is.null(cores)) {

    # detect and use half of available cores if not specified
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # run sequentially if cores is 1
  if (cores == 1L) {
    
    # run sequentially
    results_list <- lapply(

      # replications
      X = 1:reps,
      FUN = function(rep_id) {

        # one replication
        run_one_rep_study(
          rep_id        = rep_id,
          N             = N,
          T             = T,
          k             = k,
          scenarios     = scenarios,
          D_scenarios   = D_scenarios,
          A             = A,
          Psi           = Psi,
          rho_extra     = rho_extra,
          models_to_run = models_to_run,
          base_seed     = base_seed,
          ci_level      = ci_level
        )
      }
    )

    # bind and return
    return(dplyr::bind_rows(results_list))
  }

  # otherwise, run in parallel
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
      "build_riclpm_free_RI_loadings",
      "build_dpm",
      "build_dpm_free_loadings",
      "build_clpm_with_Cs",
      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_riclpm_free_RI_loadings",
      "safe_fit_dpm",
      "safe_fit_dpm_free_loadings",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "extract_lagged_parameters",
      "extract_rho_vec",
      "residualise_panel_linearC",
      "run_one_rep_study",
      "N","T","k","scenarios","D_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed","ci_level"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X  = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study(
        rep_id        = rep_id,
        N             = N,
        T             = T,
        k             = k,
        scenarios     = scenarios,
        D_scenarios   = D_scenarios,
        A             = A,
        Psi           = Psi,
        rho_extra     = rho_extra,
        models_to_run = models_to_run,
        base_seed     = base_seed,
        ci_level      = ci_level
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # bind and return
  dplyr::bind_rows(results_list)
}
