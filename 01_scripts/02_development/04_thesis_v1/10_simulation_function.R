# This is the top-level simulation runner for the simplified design:
# - one Delta trajectory per study
# - one chosen residualiser per study
# - one chosen SEM per study
#
# The function keeps the features you wanted to preserve:
# - optional one-time XGB tuning before the main replications start
# - a simple top-level interface
# - optional parallel execution
# - a progress bar through pbapply
# - a final output list that contains both the results frame and the study objects
# -------------------------------------------------------------------------------------------------

# validate study-level argument combinations before any replication starts
validate_simulation_inputs <- function(
    residualizer,
    xgb_tuning,
    tune_xgb,
    xgb_tune_args,
    residualizer_args
) {

  # match the residualiser choice
  residualizer <- match.arg(residualizer, c("none", "linear", "xgb"))

  # normalize NULL inputs to empty lists where appropriate
  if (is.null(xgb_tune_args)) {
    xgb_tune_args <- list()
  }
  if (is.null(residualizer_args)) {
    residualizer_args <- list()
  }

  # require lists
  if (!is.list(xgb_tune_args)) {
    stop("'xgb_tune_args' must be a list.")
  }
  if (!is.list(residualizer_args)) {
    stop("'residualizer_args' must be a list.")
  }

  # allowed arguments for one-time XGB tuning
  xgb_tune_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "exclude",
    "interaction_order",
    "tuning_grid",
    "cv_folds",
    "nrounds_max",
    "early_stopping_rounds",
    "nthread",
    "seed"
  )

  # allowed arguments for the linear residualiser
  linear_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "seed"
  )

  # allowed arguments for the XGB residualiser
  xgb_resid_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "nthread",
    "seed"
  )

  # no residualisation: tuning and residualiser extras must both be empty
  if (residualizer == "none") {

    if (length(xgb_tune_args) > 0) {
      stop(
        "You set residualizer = 'none', so 'xgb_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(xgb_tune_args), collapse = ", ")
      )
    }

    if (length(residualizer_args) > 0) {
      stop(
        "You set residualizer = 'none', so 'residualizer_args' must be empty. ",
        "Unused entries: ",
        paste(names(residualizer_args), collapse = ", ")
      )
    }

    if (!is.null(xgb_tuning)) {
      stop("You set residualizer = 'none', so 'xgb_tuning' must be NULL.")
    }

    if (isTRUE(tune_xgb)) {
      stop("You set residualizer = 'none', so 'tune_xgb' must be FALSE.")
    }
  }

  # linear residualisation: no XGB settings are allowed
  if (residualizer == "linear") {

    if (length(xgb_tune_args) > 0) {
      stop(
        "You set residualizer = 'linear', so 'xgb_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(xgb_tune_args), collapse = ", ")
      )
    }

    bad_linear <- setdiff(names(residualizer_args), linear_allowed)

    if (length(bad_linear) > 0) {
      stop(
        "You set residualizer = 'linear', but these are not valid linear residualiser arguments: ",
        paste(bad_linear, collapse = ", "),
        ". Allowed extra arguments are: ",
        paste(linear_allowed, collapse = ", "),
        "."
      )
    }

    if (!is.null(xgb_tuning)) {
      stop("You set residualizer = 'linear', so 'xgb_tuning' must be NULL.")
    }

    if (isTRUE(tune_xgb)) {
      stop("You set residualizer = 'linear', so 'tune_xgb' must be FALSE.")
    }
  }

  # XGB residualisation: require sensible tuning-related inputs
  if (residualizer == "xgb") {

    bad_xgb_tune <- setdiff(names(xgb_tune_args), xgb_tune_allowed)

    if (length(bad_xgb_tune) > 0) {
      stop(
        "You set residualizer = 'xgb', but these are not valid XGB tuning arguments: ",
        paste(bad_xgb_tune, collapse = ", "),
        ". Allowed tuning arguments are: ",
        paste(xgb_tune_allowed, collapse = ", "),
        "."
      )
    }

    bad_xgb_resid <- setdiff(names(residualizer_args), xgb_resid_allowed)

    if (length(bad_xgb_resid) > 0) {
      stop(
        "You set residualizer = 'xgb', but these are not valid XGB residualiser arguments: ",
        paste(bad_xgb_resid, collapse = ", "),
        ". Allowed residualiser arguments are: ",
        paste(xgb_resid_allowed, collapse = ", "),
        "."
      )
    }

    if (is.null(xgb_tuning) && !isTRUE(tune_xgb)) {
      stop(
        "You set residualizer = 'xgb'. Therefore either provide 'xgb_tuning' ",
        "or set 'tune_xgb = TRUE'."
      )
    }
  }

  invisible(TRUE)
}


run_simulation_study <- function(
    reps,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    residualizer = c("none", "linear", "xgb"),
    sem_model = c("clpm", "riclpm", "dpm"),
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE,
    bootstrap_B = 50,
    bootstrap_seed = NULL,
    cores = NULL,
    base_seed = 1234,
    xgb_tuning = NULL,
    tune_xgb = TRUE,
    xgb_tune_args = list(),
    residualizer_args = list()
) {

  # match choices
  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # validate argument combinations before doing any work
  validate_simulation_inputs(
    residualizer = residualizer,
    xgb_tuning = xgb_tuning,
    tune_xgb = tune_xgb,
    xgb_tune_args = xgb_tune_args,
    residualizer_args = residualizer_args
  )

  # choose a default number of cores
  if (is.null(cores)) {
    cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  }

  # coerce to a single positive integer
  cores <- as.integer(cores[1])

  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be a positive integer.")
  }

  # one-time XGB tuning, only if the chosen residualiser actually needs it
  if (residualizer == "xgb" && is.null(xgb_tuning) && isTRUE(tune_xgb)) {

    message("Stage 1/2: tuning XGBoost")

    # simulate one pilot data set under the same data-generating mechanism
    sim_pilot <- simulate_panel_data(
      N = N,
      T = T,
      Phi = Phi,
      Delta_list = Delta_list,
      Omega11 = Omega11,
      Sigma = Sigma,
      seed = base_seed + 1L
    )

    # newer simulator may return a list with a data component
    df_pilot <- if (is.list(sim_pilot) && !is.data.frame(sim_pilot) && !is.null(sim_pilot$data)) {
      sim_pilot$data
    } else {
      sim_pilot
    }

    # tune once and then reuse the same tuning object in every replication
    xgb_tuning <- do.call(
      tune_residualise_panel_xgb,
      c(
        list(
          df = df_pilot,
          k = k,
          exclude = exclude,
          interaction_order = confounder_order
        ),
        xgb_tune_args
      )
    )
  }

  # if XGB was requested but tuning is still missing, stop clearly
  if (residualizer == "xgb" && is.null(xgb_tuning)) {
    stop("XGB residualisation requested, but no tuning object is available.")
  }

  # create deterministic replication seeds
  rep_seeds <- base_seed + seq_len(reps)

  # run one replication
  run_one <- function(r) {
    run_one_replication(
      R = r,
      N = N,
      T = T,
      k = k,
      Phi = Phi,
      Sigma = Sigma,
      Omega11 = Omega11,
      Delta_list = Delta_list,
      residualizer = residualizer,
      sem_model = sem_model,
      confounder_order = confounder_order,
      exclude = exclude,
      free_loadings = free_loadings,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = if (is.null(bootstrap_seed)) NULL else bootstrap_seed + r,
      seed = rep_seeds[r],
      xgb_tuning = xgb_tuning,
      residualizer_args = residualizer_args
    )
  }

  # run sequentially or in parallel
  if (cores == 1L) {

    results_list <- pbapply::pblapply(
      X = seq_len(reps),
      FUN = run_one
    )

  } else {

    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(lavaan)
      library(pbapply)
      library(xgboost)
      NULL
    })

    parallel::clusterExport(
      cl = cl,
      varlist = ls(envir = .GlobalEnv),
      envir = .GlobalEnv
    )

    results_list <- pbapply::pblapply(
      X = seq_len(reps),
      FUN = run_one,
      cl = cl
    )
  }

  # stack all replication outputs
  results <- do.call(rbind, results_list)

  # return both the final results and the study objects
  list(
    results = results,
    study = list(
      reps = reps,
      N = N,
      T = T,
      k = k,
      Phi = Phi,
      Sigma = Sigma,
      Omega11 = Omega11,
      Delta_list = Delta_list,
      residualizer = residualizer,
      sem_model = sem_model,
      confounder_order = confounder_order,
      exclude = exclude,
      free_loadings = free_loadings,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      cores = cores,
      base_seed = base_seed,
      xgb_tuning = xgb_tuning,
      tune_xgb = tune_xgb,
      xgb_tune_args = xgb_tune_args,
      residualizer_args = residualizer_args
    )
  )
}