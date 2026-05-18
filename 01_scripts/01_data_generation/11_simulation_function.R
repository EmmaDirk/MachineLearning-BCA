# =================================================================================================
#
# This script contains the shared simulation engine used by the top-level runner.
#
# The efficient model-set interface lets you define multiple analysis pipelines
# and then runs them in the computationally cheapest order:
# - one simulated data set per replication
# - one prepared stage-1 data set per unique residualiser recipe
# - all compatible SEMs fit on that same prepared data
# - bootstrap draws shared the same way
#
# Single-model use is supported through a wrapper around the same shared engine.
#
# Thread and worker counts are controlled by the top-level simulation runner. This script only
# receives the requested number of worker processes through the cores argument.
# =================================================================================================

# ---- validation helpers --------------------------------------------------------------------------

# Validate one model specification before any replication starts.

validate_simulation_inputs <- function(
    residualizer,
    xgb_tuning,
    enet_tuning,
    tune_xgb,
    tune_enet,
    xgb_tune_args,
    enet_tune_args,
    residualizer_args
) {

  # Match the residualiser choice.

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))

  # Convert NULL argument lists to empty lists.

  if (is.null(xgb_tune_args)) {
    xgb_tune_args <- list()
  }

  if (is.null(enet_tune_args)) {
    enet_tune_args <- list()
  }

  if (is.null(residualizer_args)) {
    residualizer_args <- list()
  }

  # Require all optional argument containers to be lists.

  if (!is.list(xgb_tune_args)) {
    stop("'xgb_tune_args' must be a list.")
  }

  if (!is.list(enet_tune_args)) {
    stop("'enet_tune_args' must be a list.")
  }

  if (!is.list(residualizer_args)) {
    stop("'residualizer_args' must be a list.")
  }

  # ---- allowed argument names --------------------------------------------------------------------

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
    "early_tune_wave",
    "late_tune_wave",
    "late_start_wave",
    "nthread",
    "seed"
  )

  enet_tune_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "exclude",
    "interaction_order",
    "alpha_grid",
    "cv_folds",
    "seed"
  )

  linear_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "seed"
  )

  xgb_resid_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "use_oof",
    "late_start_wave",
    "nthread",
    "seed"
  )

  enet_resid_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "seed"
  )

  # ---- no residualisation ------------------------------------------------------------------------

  if (residualizer == "none") {

    if (length(xgb_tune_args) > 0) {
      stop(
        "You set residualizer = 'none', so 'xgb_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(xgb_tune_args), collapse = ", ")
      )
    }

    if (length(enet_tune_args) > 0) {
      stop(
        "You set residualizer = 'none', so 'enet_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(enet_tune_args), collapse = ", ")
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

    if (!is.null(enet_tuning)) {
      stop("You set residualizer = 'none', so 'enet_tuning' must be NULL.")
    }

    if (isTRUE(tune_xgb)) {
      stop("You set residualizer = 'none', so 'tune_xgb' must be FALSE.")
    }

    if (isTRUE(tune_enet)) {
      stop("You set residualizer = 'none', so 'tune_enet' must be FALSE.")
    }
  }

  # ---- linear residualisation --------------------------------------------------------------------

  if (residualizer == "linear") {

    if (length(xgb_tune_args) > 0) {
      stop(
        "You set residualizer = 'linear', so 'xgb_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(xgb_tune_args), collapse = ", ")
      )
    }

    if (length(enet_tune_args) > 0) {
      stop(
        "You set residualizer = 'linear', so 'enet_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(enet_tune_args), collapse = ", ")
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

    if (!is.null(enet_tuning)) {
      stop("You set residualizer = 'linear', so 'enet_tuning' must be NULL.")
    }

    if (isTRUE(tune_xgb)) {
      stop("You set residualizer = 'linear', so 'tune_xgb' must be FALSE.")
    }

    if (isTRUE(tune_enet)) {
      stop("You set residualizer = 'linear', so 'tune_enet' must be FALSE.")
    }
  }

  # ---- XGB residualisation -----------------------------------------------------------------------

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

    if (!is.null(enet_tuning)) {
      stop("You set residualizer = 'xgb', so 'enet_tuning' must be NULL.")
    }

    if (length(enet_tune_args) > 0) {
      stop("You set residualizer = 'xgb', so 'enet_tune_args' must be empty.")
    }

    if (isTRUE(tune_enet)) {
      stop("You set residualizer = 'xgb', so 'tune_enet' must be FALSE.")
    }
  }

  # ---- elastic net residualisation ---------------------------------------------------------------

  if (residualizer == "enet") {

    if (length(xgb_tune_args) > 0) {
      stop(
        "You set residualizer = 'enet', so 'xgb_tune_args' must be empty. ",
        "Unused entries: ",
        paste(names(xgb_tune_args), collapse = ", ")
      )
    }

    bad_enet_tune <- setdiff(names(enet_tune_args), enet_tune_allowed)

    if (length(bad_enet_tune) > 0) {
      stop(
        "You set residualizer = 'enet', but these are not valid Elastic Net tuning arguments: ",
        paste(bad_enet_tune, collapse = ", "),
        ". Allowed tuning arguments are: ",
        paste(enet_tune_allowed, collapse = ", "),
        "."
      )
    }

    bad_enet_resid <- setdiff(names(residualizer_args), enet_resid_allowed)

    if (length(bad_enet_resid) > 0) {
      stop(
        "You set residualizer = 'enet', but these are not valid Elastic Net residualiser arguments: ",
        paste(bad_enet_resid, collapse = ", "),
        ". Allowed residualiser arguments are: ",
        paste(enet_resid_allowed, collapse = ", "),
        "."
      )
    }

    if (is.null(enet_tuning) && !isTRUE(tune_enet)) {
      stop(
        "You set residualizer = 'enet'. Therefore either provide 'enet_tuning' ",
        "or set 'tune_enet = TRUE'."
      )
    }

    if (!is.null(xgb_tuning)) {
      stop("You set residualizer = 'enet', so 'xgb_tuning' must be NULL.")
    }

    if (isTRUE(tune_xgb)) {
      stop("You set residualizer = 'enet', so 'tune_xgb' must be FALSE.")
    }
  }

  invisible(TRUE)
}


# Validate every model specification before any replication starts.

validate_model_spec_list <- function(model_specs) {

  for (spec in model_specs) {
    validate_simulation_inputs(
      residualizer = spec$residualizer,
      xgb_tuning = spec$xgb_tuning,
      enet_tuning = spec$enet_tuning,
      tune_xgb = spec$tune_xgb,
      tune_enet = spec$tune_enet,
      xgb_tune_args = spec$xgb_tune_args,
      enet_tune_args = spec$enet_tune_args,
      residualizer_args = spec$residualizer_args
    )
  }

  invisible(TRUE)
}


# ---- worker script resolution -------------------------------------------------------------------

# Define the default script order used when workers source the simulation code.

default_worker_script_files <- function() {

  c(
    "01_delta_sampler.R",
    "02_delta_trajectory.R",
    "03_simulate_panel_data.R",
    "04_lavaan_model_string_builder.R",
    "05_residualisers.R",
    "06_model_fitters.R",
    "07_bootstrap_helpers.R",
    "08_fit_stat_extractors.R",
    "09_model_set_helpers.R",
    "10_one_replication_wrapper.R",
    "11_simulation_function.R"
  )
}


# Resolve the worker script file set once, so master and workers use the same code path.

resolve_worker_script_files <- function(script_dir, script_files = NULL) {

  if (is.null(script_dir)) {
    return(NULL)
  }

  if (is.null(script_files)) {
    script_files <- default_worker_script_files()
  }

  script_files <- as.character(script_files)

  missing_files <- script_files[!file.exists(file.path(script_dir, script_files))]

  if (length(missing_files) > 0L) {
    stop(
      "The following worker script files were requested but not found: ",
      paste(missing_files, collapse = ", ")
    )
  }

  script_files
}


# ---- worker helper export -----------------------------------------------------------------------

# Collect the simulation helper functions that should be available on every worker.

simulation_helper_function_names <- function(helper_env = environment(run_simulation_model_set)) {

  helper_names <- ls(envir = helper_env, all.names = FALSE)

  helper_names[vapply(helper_names, function(nm) {
    is.function(get(nm, envir = helper_env, inherits = FALSE))
  }, logical(1))]
}


# Export helper functions directly when worker-side sourcing is not available.

export_simulation_helpers_to_cluster <- function(
    cl,
    helper_env = environment(run_simulation_model_set),
    helper_names = NULL
) {

  if (is.null(helper_names)) {
    helper_names <- simulation_helper_function_names(helper_env = helper_env)
  }

  if (length(helper_names) == 0) {
    return(invisible(FALSE))
  }

  parallel::clusterExport(
    cl = cl,
    varlist = helper_names,
    envir = helper_env
  )

  invisible(TRUE)
}


# Export only the objects required by the current run.

export_run_context_to_cluster <- function(
    cl,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    model_specs,
    stage1_groups,
    burn_in,
    bootstrap_seed,
    rep_seeds
) {

  parallel::clusterExport(
    cl = cl,
    varlist = c(
      "N",
      "T",
      "k",
      "Phi",
      "Sigma",
      "Omega11",
      "Delta_list",
      "model_specs",
      "stage1_groups",
      "burn_in",
      "bootstrap_seed",
      "rep_seeds"
    ),
    envir = environment()
  )

  invisible(TRUE)
}


# ---- worker entry point -------------------------------------------------------------------------

# Run one replication from the context exported to the worker.

run_one_replication_from_exported_context <- function(r) {

  r <- as.integer(r[1])

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


# ---- cluster setup ------------------------------------------------------------------------------

# Initialise a PSOCK cluster without exporting the whole global environment.

make_simulation_cluster <- function(
    cores,
    script_dir = NULL,
    script_files = NULL,
    helper_env = environment(run_simulation_model_set)
) {

  cl <- parallel::makeCluster(cores)

  # Load required packages on each worker. Thread limits are set by the top-level runner.

  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(pbapply)
    library(xgboost)
    library(glmnet)

    NULL
  })

  resolved_script_files <- resolve_worker_script_files(
    script_dir = script_dir,
    script_files = script_files
  )

  # Source the exact same worker script set when available.

  if (!is.null(script_dir) && !is.null(resolved_script_files)) {

    parallel::clusterExport(
      cl = cl,
      varlist = c("script_dir", "resolved_script_files"),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      for (f in resolved_script_files) {
        source(file.path(script_dir, f))
      }

      NULL
    })

  } else {

    # Export helper functions only when no script directory is available.

    export_simulation_helpers_to_cluster(
      cl = cl,
      helper_env = helper_env
    )
  }

  cl
}


# ---- XGB tuning ---------------------------------------------------------------------------------

# Tune XGB once per unique XGB recipe and assign the result back to every matching model.

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

  if (length(xgb_specs) == 0) {
    return(specs)
  }

  needs_internal_tuning <- any(vapply(
    xgb_specs,
    function(x) is.null(x$xgb_tuning) && isTRUE(x$tune_xgb),
    logical(1)
  ))

  df_pilot <- NULL

  if (needs_internal_tuning) {

    message("Stage 1/2: tuning XGBoost")

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

  tuning_group_ids <- sort(unique(vapply(
    xgb_specs,
    function(x) x$xgb_tuning_group_id,
    integer(1)
  )))

  for (gid in tuning_group_ids) {

    idx <- which(vapply(
      specs,
      function(x) identical(x$xgb_tuning_group_id, gid),
      logical(1)
    ))

    proto <- specs[[idx[1]]]
    tuning_obj <- proto$xgb_tuning

    if (is.null(tuning_obj) && isTRUE(proto$tune_xgb)) {

      tuning_obj <- do.call(
        tune_residualise_panel_xgb,
        c(
          list(
            df = df_pilot,
            k = k,
            exclude = proto$residualizer_exclude,
            interaction_order = proto$residualizer_c_order
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


# ---- elastic net tuning -------------------------------------------------------------------------

# Tune Elastic Net once per unique Elastic Net recipe and assign the result back to each model.

resolve_enet_tuning_objects <- function(
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

  specs <- assign_enet_tuning_group_ids(model_specs)

  enet_specs <- Filter(function(x) x$residualizer == "enet", specs)

  if (length(enet_specs) == 0) {
    return(specs)
  }

  needs_internal_tuning <- any(vapply(
    enet_specs,
    function(x) is.null(x$enet_tuning) && isTRUE(x$tune_enet),
    logical(1)
  ))

  df_pilot <- NULL

  if (needs_internal_tuning) {

    message("Stage 1.5/2: tuning Elastic Net")

    df_pilot <- simulate_panel_data(
      N = N,
      T = T,
      Phi = Phi,
      Delta_list = Delta_list,
      Omega11 = Omega11,
      Sigma = Sigma,
      burn_in = burn_in,
      seed = base_seed + 2L
    )
  }

  tuning_group_ids <- sort(unique(vapply(
    enet_specs,
    function(x) x$enet_tuning_group_id,
    integer(1)
  )))

  for (gid in tuning_group_ids) {

    idx <- which(vapply(
      specs,
      function(x) identical(x$enet_tuning_group_id, gid),
      logical(1)
    ))

    proto <- specs[[idx[1]]]
    tuning_obj <- proto$enet_tuning

    if (is.null(tuning_obj) && isTRUE(proto$tune_enet)) {

      tuning_obj <- do.call(
        tune_residualise_panel_enet,
        c(
          list(
            df = df_pilot,
            k = k,
            exclude = proto$residualizer_exclude,
            interaction_order = proto$residualizer_c_order
          ),
          proto$enet_tune_args
        )
      )
    }

    for (i in idx) {
      specs[[i]]$enet_tuning <- tuning_obj
    }
  }

  specs
}


# ---- shared simulation engine -------------------------------------------------------------------

# Run the efficient multi-model simulation study.

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

  # Check the burn-in length.

  if (!is.numeric(burn_in) || length(burn_in) != 1 || is.na(burn_in) ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("'burn_in' must be a single non-negative integer.")
  }

  burn_in <- as.integer(burn_in)

  # Validate the full model set before running any replication.

  model_specs <- normalize_model_spec_list(model_specs)
  validate_model_spec_list(model_specs)

  # Require the top-level runner to define the number of workers.

  cores <- as.integer(cores[1])

  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be supplied by the top-level runner as a positive integer.")
  }

  # Resolve tuning objects before building the stage-1 groups.

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

  model_specs <- resolve_enet_tuning_objects(
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

  # Build the shared stage-1 execution plan.

  model_specs <- assign_stage1_group_ids(model_specs)
  stage1_groups <- build_stage1_groups(model_specs)

  # Create deterministic replication seeds.

  rep_seeds <- base_seed + seq_len(reps)

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

  # Run sequentially or in parallel.

  if (cores == 1L) {

    results_list <- pbapply::pblapply(
      X = seq_len(reps),
      FUN = run_one
    )

  } else {

    script_dir_local <- NULL

    if (exists("script_dir", inherits = TRUE)) {
      script_dir_local <- get("script_dir", inherits = TRUE)
    }

    script_files_local <- NULL

    if (exists("simulation_worker_script_files", inherits = TRUE)) {
      script_files_local <- get("simulation_worker_script_files", inherits = TRUE)
    } else if (exists("simulation_script_files", inherits = TRUE)) {
      script_files_local <- get("simulation_script_files", inherits = TRUE)
    }

    cl <- make_simulation_cluster(
      cores = cores,
      script_dir = script_dir_local,
      script_files = script_files_local,
      helper_env = environment(run_simulation_model_set)
    )

    on.exit(parallel::stopCluster(cl), add = TRUE)

    export_run_context_to_cluster(
      cl = cl,
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
      bootstrap_seed = bootstrap_seed,
      rep_seeds = rep_seeds
    )

    results_list <- pbapply::pblapply(
      X = seq_len(reps),
      FUN = run_one_replication_from_exported_context,
      cl = cl
    )
  }

  # Stack all replication outputs and split them by model name.

  results <- do.call(rbind, results_list)
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


# ---- single-model wrapper -----------------------------------------------------------------------

# Run a single model through the shared model-set engine.

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
    residualizer = c("none", "linear", "xgb", "enet"),
    sem_model = c("clpm", "riclpm", "dpm"),
    sem_c_order = 0,
    sem_exclude = NULL,
    residualizer_c_order = 0,
    residualizer_exclude = NULL,
    free_loadings = FALSE,
    bootstrap_B = 50,
    bootstrap_seed = NULL,
    cores = NULL,
    base_seed = 1234,
    xgb_tuning = NULL,
    enet_tuning = NULL,
    tune_xgb = TRUE,
    tune_enet = TRUE,
    xgb_tune_args = list(),
    enet_tune_args = list(),
    residualizer_args = list()
) {

  # Match the requested model choices.

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # Create one model specification and run it through the shared engine.

  spec <- make_model_spec(
    name = "model_1",
    residualizer = residualizer,
    sem_model = sem_model,
    sem_c_order = sem_c_order,
    sem_exclude = sem_exclude,
    residualizer_c_order = residualizer_c_order,
    residualizer_exclude = residualizer_exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    xgb_tuning = xgb_tuning,
    enet_tuning = enet_tuning,
    tune_xgb = tune_xgb,
    tune_enet = tune_enet,
    xgb_tune_args = xgb_tune_args,
    enet_tune_args = enet_tune_args,
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