# This is the top-level simulation runner.
#
# The efficient model-set interface lets you define multiple analysis pipelines
# and then runs them in the computationally cheapest order:
# - one simulated data set per replication
# - one prepared stage-1 data set per unique residualiser recipe
# - all compatible SEMs fit on that same prepared data
# - bootstrap draws shared the same way
#
# Single-model use is still supported, but it is routed through the same shared engine.
# -------------------------------------------------------------------------------------------------

# validate one single model-spec argument combination before any replication starts
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


# validate every model specification before any replication starts
validate_model_spec_list <- function(model_specs) {

  for (spec in model_specs) {
    validate_simulation_inputs(
      residualizer = spec$residualizer,
      xgb_tuning = spec$xgb_tuning,
      tune_xgb = spec$tune_xgb,
      xgb_tune_args = spec$xgb_tune_args,
      residualizer_args = spec$residualizer_args
    )
  }

  invisible(TRUE)
}


# choose a sensible number of worker cores
resolve_core_count <- function(cores) {

  if (is.null(cores)) {
    cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  }

  cores <- as.integer(cores[1])

  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be a positive integer.")
  }

  cores
}


# tune XGB only once per unique XGB tuning recipe and assign the result back to every model spec
resolve_xgb_tuning_objects <- function(
    model_specs,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    burn_in,
    base_seed
) {

  specs <- assign_xgb_tuning_group_ids(model_specs)

  xgb_specs <- Filter(function(x) x$residualizer == "xgb", specs)

  # if no XGB model is present, nothing to do
  if (length(xgb_specs) == 0) {
    return(specs)
  }

  # do we need to run internal tuning at all?
  needs_internal_tuning <- any(vapply(
    xgb_specs,
    function(x) is.null(x$xgb_tuning) && isTRUE(x$tune_xgb),
    logical(1)
  ))

  df_pilot <- NULL

  if (needs_internal_tuning) {

    message("Stage 1/2: tuning XGBoost")

    # simulate one pilot data set under the same data-generating mechanism
    df_pilot <- simulate_panel_data(
      N = N,
      T = T,
      Phi = Phi,
      Delta_list = Delta_list,
      Omega11 = Omega11,
      Sigma = Sigma,
      burn_in = burn_in,
      seed = base_seed + 1L
    )
  }

  # resolve one tuning object per unique tuning group
  tuning_group_ids <- sort(unique(vapply(xgb_specs, function(x) x$xgb_tuning_group_id, integer(1))))

  for (gid in tuning_group_ids) {

    idx <- which(vapply(specs, function(x) identical(x$xgb_tuning_group_id, gid), logical(1)))
    proto <- specs[[idx[1]]]

    tuning_obj <- proto$xgb_tuning

    if (is.null(tuning_obj) && isTRUE(proto$tune_xgb)) {

      tuning_obj <- do.call(
        tune_residualise_panel_xgb,
        c(
          list(
            df = df_pilot,
            k = k,
            exclude = proto$exclude,
            interaction_order = proto$confounder_order
          ),
          proto$xgb_tune_args
        )
      )
    }

    for (i in idx) {
      specs[[i]]$xgb_tuning <- tuning_obj
    }
  }

  specs
}


# run the efficient multi-model simulation study
run_simulation_model_set <- function(
    reps,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    model_specs,
    burn_in = 0L,
    bootstrap_seed = NULL,
    cores = NULL,
    base_seed = 1234
) {

  # check burn-in
  if (!is.numeric(burn_in) || length(burn_in) != 1 || is.na(burn_in) ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("'burn_in' must be a single non-negative integer.")
  }
  burn_in <- as.integer(burn_in)

  # normalize and validate the whole model set
  model_specs <- normalize_model_spec_list(model_specs)
  validate_model_spec_list(model_specs)

  # choose a default number of cores
  cores <- resolve_core_count(cores)

  # resolve XGB tuning objects before we build the stage-1 groups
  model_specs <- resolve_xgb_tuning_objects(
    model_specs = model_specs,
    N = N,
    T = T,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    burn_in = burn_in,
    base_seed = base_seed
  )

  # build the final shared stage-1 execution plan
  model_specs <- assign_stage1_group_ids(model_specs)
  stage1_groups <- build_stage1_groups(model_specs)

  # create deterministic replication seeds
  rep_seeds <- base_seed + seq_len(reps)

  # run one replication
  run_one <- function(r) {
    run_one_replication_model_set(
      R = r,
      N = N,
      T = T,
      k = k,
      Phi = Phi,
      Sigma = Sigma,
      Omega11 = Omega11,
      Delta_list = Delta_list,
      model_specs = model_specs,
      stage1_groups = stage1_groups,
      burn_in = burn_in,
      bootstrap_seed = if (is.null(bootstrap_seed)) NULL else bootstrap_seed + r,
      seed = rep_seeds[r]
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

  # also return one data frame per originally requested model
  results_by_model <- split_results_by_model_name(results, model_specs)

  list(
    results = results,
    results_by_model = results_by_model,
    study = list(
      reps = reps,
      N = N,
      T = T,
      k = k,
      Phi = Phi,
      Sigma = Sigma,
      Omega11 = Omega11,
      Delta_list = Delta_list,
      burn_in = burn_in,
      bootstrap_seed = bootstrap_seed,
      cores = cores,
      base_seed = base_seed,
      model_specs = model_specs,
      stage1_groups = stage1_groups
    )
  )
}


# convenience wrapper for running a single model through the shared model-set engine
run_simulation_study <- function(
    reps,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    burn_in = 0L,
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

  # create one model specification and hand it to the shared engine
  spec <- make_model_spec(
    name = "model_1",
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    xgb_tuning = xgb_tuning,
    tune_xgb = tune_xgb,
    xgb_tune_args = xgb_tune_args,
    residualizer_args = residualizer_args
  )

  run_simulation_model_set(
    reps = reps,
    N = N,
    T = T,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    model_specs = list(spec),
    burn_in = burn_in,
    bootstrap_seed = bootstrap_seed,
    cores = cores,
    base_seed = base_seed
  )
}