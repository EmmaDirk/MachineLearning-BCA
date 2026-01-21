# this script runs the wrapper function that executes one replication of the simulation study
# adding a progress bar and parallelization
# ------------------------------------------------------------------------------------------------

run_simulation_study <- function(
  reps,                                                                    # number of replications
  N,                                                                       # sample size
  T,                                                                       # number of waves
  k,                                                                       # number of linear confounders
  scenarios,                                                               # e.g., c("constant","stepwise")
  B_scenarios,                                                             # user-supplied B matrices per scenario
  A,                                                                       # 2×2 AR + cross-lag matrix
  Psi,                                                                     # k×k confounder covariance
  rho_extra,                                                               # extra covariance added to observations
  models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca","xgb")
  cores = NULL,                                                            # default is detectCores()/2
  base_seed = 1234,                                                        # master seed for reproducible reps
  ci_level = 0.95                                                          # CI level for extracted parameters
) {

  # if cores is not specified, detect and use half of available cores
  if (is.null(cores)) {
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # XGB one-time tuning for this study
  if ("xgb" %in% models_to_run) {

    # pick a pilot scenario for tuning
    scen_tune <- scenarios[1]

    # set seed for the pilot
    set.seed(base_seed + 1)

    # check B_scenarios is a named list and contains the pilot scenario
    if (!is.list(B_scenarios) || is.null(names(B_scenarios)) || !scen_tune %in% names(B_scenarios))
      stop("For xgb tuning, B_scenarios must be a named list containing the first scenario name.")

    # take the pilot B trajectory
    B_list_pilot <- B_scenarios[[scen_tune]]

    # check pilot trajectory is a list of length T
    if (!is.list(B_list_pilot))
      stop("For xgb tuning, B_scenarios[['", scen_tune, "']] must be a list of length T.")
    if (length(B_list_pilot) != T)
      stop("For xgb tuning, B_scenarios[['", scen_tune, "']] has length ", length(B_list_pilot), " but T = ", T, ".")

    # simulate pilot data
    df_pilot <- simulate_panel_data_int(
      N         = N,
      T         = T,
      A         = A,
      B_list    = B_list_pilot,
      Psi       = Psi,
      rho_extra = rho_extra
    )

    # tune and store for residualise_panel_xgb
    XGB_TUNED <- tune_xgb_once(df_pilot, k)
    assign("XGB_TUNED", XGB_TUNED, envir = .GlobalEnv)

  } else {

    # set tuning object to NULL if xgb is not used
    XGB_TUNED <- NULL
    assign("XGB_TUNED", XGB_TUNED, envir = .GlobalEnv)
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
          base_seed     = base_seed,
          ci_level      = ci_level
        )
      }
    )

    return(dplyr::bind_rows(results_list))
  }

  # make the cluster
  cl <- parallel::makeCluster(cores)

  # load packages on each worker
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(mvtnorm)
    library(xgboost)
    NULL
  })

  # export functions and objects needed on workers
  parallel::clusterExport(
    cl,
    c(
      "simulate_panel_data_int",

      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",

      "residualise_panel_linearC",
      "residualise_panel_xgb",

      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "safe_fit_clpm_xgb",

      "extract_lagged_parameters",
      "extract_rho_vec",

      "run_1_rep_study",

      "N","T","k","scenarios","B_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed","ci_level",

      "XGB_TUNED"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X  = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_1_rep_study(
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
        base_seed     = base_seed,
        ci_level      = ci_level
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # return results
  dplyr::bind_rows(results_list)
}
