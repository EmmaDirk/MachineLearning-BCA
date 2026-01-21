# 10_simulation_function: # this script runs the wrapper function that executes one replication of the simulation study
# adding a progress bar and parallelization
#
# This function runs the full simulation study by repeatedly calling run_one_rep_study() across replications,
# either sequentially (single core) or in parallel (multiple cores).
#
# 1) It determines how many CPU cores to use:
# - if cores is NULL, it defaults to detectCores()/2 (at least 1)
# - if cores == 1, it runs everything sequentially via lapply()
#
# 2) It handles one-time XGB tuning (only if requested):
# - if "xgb" is in models_to_run and xgb_tuned is NULL and tune_xgb is TRUE, it tunes once on a pilot dataset
# - the pilot dataset is simulated from the first scenario using base_seed + 1
# - the tuned object is stored in xgb_tuned and passed into run_one_rep_study()
#
# 3) If running sequentially (cores == 1):
# - for rep_id = 1..reps, it calls run_one_rep_study() with the provided settings (including xgb_tuned)
# - it row-binds the per-replication outputs into one combined results data frame
#
# 4) If running in parallel (cores > 1):
# - it creates a PSOCK cluster with the requested number of workers
# - it initializes reproducible random number streams on the cluster (clusterSetRNGStream with base_seed)
# - it loads required packages on each worker (lavaan, mvtnorm, xgboost)
# - it exports all functions and objects needed by the workers (including D_scenarios, model settings, and xgb_tuned)
#
# 5) It distributes replications across workers:
# - it uses pbapply::pblapply() to run rep_id = 1..reps in parallel with a progress bar
# - each worker calls run_one_rep_study(), which handles scenario checks, data simulation, model fitting, and extraction
#
# 6) It finalizes and returns results:
# - it stops the cluster to free resources
# - it row-binds all replication outputs into a single results data frame and returns it
# ------------------------------------------------------------------------------------------------------------

run_simulation_study <- function(
  reps,                                                                    # number of replications
  N,                                                                       # sample size
  T,                                                                       # number of waves
  k,                                                                       # number of linear confounders
  scenarios,                                                               # e.g., c("constant","stepwise")
  D_scenarios,                                                             # user-supplied delta matrices per scenario
  A,                                                                       # 2×2 AR + cross-lag matrix
  Psi,                                                                     # k×k confounder covariance
  rho_extra,                                                               # extra covariance added to observations
  models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca","xgb")
  cores = NULL,                                                            # default is detectCores()/2
  base_seed = 1234,                                                        # master seed for reproducible reps
  ci_level = 0.95,                                                         # CI level for extracted parameters
  xgb_tuned = NULL,                                                        # optionally provide a tuned object directly
  tune_xgb = TRUE,                                                         # if TRUE and xgb_tuned is NULL, tune once on a pilot dataset
  xgb_tune_profile = c("quick", "medium", "full"),
  xgb_fit_profile = c("fast", "balanced", "thorough"),
  xgb_tune_grid = NULL,
  xgb_tune_overrides = NULL,
  xgb_fit_overrides = NULL,
  xgb_profile = NULL
) {

  # xgb tuning settings are passed as arguments
  # xgb fitting settings are passed as arguments

  # allow legacy alias from the master script
  if (!is.null(xgb_profile)) {
    xgb_tune_profile <- xgb_profile
  }

  # match xgb profile arguments
  xgb_tune_profile <- match.arg(xgb_tune_profile)
  xgb_fit_profile  <- match.arg(xgb_fit_profile)

  # if cores is not specified, detect and use half of available cores
  if (is.null(cores)) {
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # one-time XGB tuning for this study (optional)
  if ("xgb" %in% models_to_run && is.null(xgb_tuned) && isTRUE(tune_xgb)) {

    # pick a pilot scenario for tuning
    scen_tune <- scenarios[1]

    # set seed for the pilot
    set.seed(base_seed + 1)

    # check D_scenarios is a named list and contains the pilot scenario
    if (!is.list(D_scenarios) || is.null(names(D_scenarios)) || !scen_tune %in% names(D_scenarios))
      stop("For xgb tuning, D_scenarios must be a named list containing the first scenario name.")

    # take the pilot delta trajectory
    D_list_pilot <- D_scenarios[[scen_tune]]

    # check pilot trajectory is a list of length T
    if (!is.list(D_list_pilot))
      stop("For xgb tuning, D_scenarios[['", scen_tune, "']] must be a list of length T.")
    if (length(D_list_pilot) != T)
      stop("For xgb tuning, D_scenarios[['", scen_tune, "']] has length ", length(D_list_pilot), " but T = ", T, ".")

    # simulate pilot data
    df_pilot <- simulate_panel_data_int(
      N         = N,
      T         = T,
      A         = A,
      D_list    = D_list_pilot,
      Psi       = Psi,
      rho_extra = rho_extra
    )

    # build tuning arguments
    tune_args <- list(
      df = df_pilot,
      k = k,
      tune_profile = xgb_tune_profile,
      tune_grid = xgb_tune_grid,
      seed = 1
    )

    # allow overwrite options for tuning
    if (!is.null(xgb_tune_overrides)) {
      if (!is.list(xgb_tune_overrides)) stop("xgb_tune_overrides must be a list when provided.")
      tune_args <- c(tune_args, xgb_tune_overrides)
    }

    # tune and store locally for passing into run_one_rep_study
    xgb_tuned <- do.call(tune_xgb_once, tune_args)
  }

  # run sequentially if cores is 1
  if (cores == 1L) {

    results_list <- lapply(
      X = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study(
          rep_id          = rep_id,
          N               = N,
          T               = T,
          k               = k,
          scenarios       = scenarios,
          B_scenarios     = D_scenarios,
          A               = A,
          Psi             = Psi,
          rho_extra       = rho_extra,
          models_to_run   = models_to_run,
          base_seed       = base_seed,
          ci_level        = ci_level,
          xgb_tuned       = xgb_tuned,
          xgb_fit_profile = xgb_fit_profile,
          xgb_fit_overrides = xgb_fit_overrides
        )
      }
    )

    return(dplyr::bind_rows(results_list))
  }

  # make the cluster
  cl <- parallel::makeCluster(cores)

  # initialize reproducible RNG streams on the cluster
  parallel::clusterSetRNGStream(cl, iseed = base_seed)

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
      "tune_xgb_once",

      "build_clpm",
      "build_clpm_with_Cs",
      "build_riclpm",
      "build_dpm",

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

      "n_interactions_from_k",
      "run_one_rep_study",

      "N","T","k","scenarios","D_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed","ci_level",

      "xgb_tuned",
      "xgb_fit_profile",
      "xgb_fit_overrides"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X  = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study(
        rep_id          = rep_id,
        N               = N,
        T               = T,
        k               = k,
        scenarios       = scenarios,
        B_scenarios     = D_scenarios,
        A               = A,
        Psi             = Psi,
        rho_extra       = rho_extra,
        models_to_run   = models_to_run,
        base_seed       = base_seed,
        ci_level        = ci_level,
        xgb_tuned       = xgb_tuned,
        xgb_fit_profile = xgb_fit_profile,
        xgb_fit_overrides = xgb_fit_overrides
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # return results
  dplyr::bind_rows(results_list)
}
