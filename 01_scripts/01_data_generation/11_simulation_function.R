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
    enet_tuning,
    tune_xgb,
    tune_enet,
    xgb_tune_args,
    enet_tune_args,
    residualizer_args
) {

  # match the residualiser choice
  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))

  # normalize NULL inputs to empty lists where appropriate
  if (is.null(xgb_tune_args)) {
    xgb_tune_args <- list()
  }
  if (is.null(enet_tune_args)) {
    enet_tune_args <- list()
  }
  if (is.null(residualizer_args)) {
    residualizer_args <- list()
  }

  # require lists
  if (!is.list(xgb_tune_args)) {
    stop("'xgb_tune_args' must be a list.")
  }
  if (!is.list(enet_tune_args)) {
    stop("'enet_tune_args' must be a list.")
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
    "early_tune_wave",
    "late_tune_wave",
    "late_start_wave",
    "nthread",
    "seed"
  )

  # allowed arguments for Elastic Net one-time tuning
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
    "use_oof",
    "late_start_wave",
    "nthread",
    "seed"
  )

  # allowed arguments for the Elastic Net residualiser
  enet_resid_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
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

  # linear residualisation: no XGB settings are allowed
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

  # Elastic Net residualisation: require sensible tuning-related inputs
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


# validate every model specification before any replication starts
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


# apply conservative thread caps for threaded math libraries inside each R process
# This is critical when we parallelise at the replication level, because otherwise each
# worker process may itself try to use many BLAS / OpenMP threads.
set_single_thread_math_env <- function() {

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OMP_DYNAMIC = "FALSE",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    MKL_DYNAMIC = "FALSE",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  invisible(TRUE)
}


# best-effort helper to also cap BLAS / OpenMP threads from inside R when RhpcBLASctl
# is available. The simulation still works when the package is not installed.
cap_threads_runtime <- function() {

  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)

  invisible(TRUE)
}


# default script file order for worker-side sourcing
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


# resolve the exact worker script file set once so the master and workers stay on the same code path
resolve_worker_script_files <- function(script_dir, script_files = NULL) {

  if (is.null(script_dir)) {
    return(NULL)
  }

  # if the caller already passed an explicit file vector, use exactly that vector
  if (!is.null(script_files)) {
    script_files <- as.character(script_files)

    missing_files <- script_files[!file.exists(file.path(script_dir, script_files))]
    if (length(missing_files) > 0) {
      stop(
        "The following worker script files were requested but not found: ",
        paste(missing_files, collapse = ", ")
      )
    }

    return(script_files)
  }

  # otherwise, choose a canonical set with backward-compatible fallbacks for renamed files
  candidate_groups <- list(
    c("01_delta_sampler.R"),
    c("02_delta_trajectory.R"),
    c("03_simulate_panel_data.R"),
    c("04_lavaan_model_string_builder.R"),
    c("05_residualisers.R"),
    c("06_model_fitters.R"),
    c("07_bootstrap_helpers.R", "07_bootstrap_helpers_fixed.R", "07_bootstrap_helpers_nobc.R"),
    c("08_fit_stat_extractors.R"),
    c("09_model_set_helpers.R", "09_model_set_helpers_nobc.R"),
    c("10_one_replication_wrapper.R", "10_one_replication_wrapper_nobc.R"),
    c("11_simulation_function.R", "11_simulation_function_nobc.R")
  )

  chosen <- vapply(candidate_groups, function(candidates) {
    hit <- candidates[file.exists(file.path(script_dir, candidates))]

    if (length(hit) == 0) {
      stop(
        "Could not resolve a required worker script in ", script_dir,
        ". Tried: ", paste(candidates, collapse = ", ")
      )
    }

    hit[1]
  }, character(1))

  as.character(chosen)
}


# collect the simulation helper functions that should be available on every worker
simulation_helper_function_names <- function(helper_env = environment(run_simulation_model_set)) {

  helper_names <- ls(envir = helper_env, all.names = FALSE)

  helper_names[vapply(helper_names, function(nm) {
    is.function(get(nm, envir = helper_env, inherits = FALSE))
  }, logical(1))]
}


# export the helper functions directly when worker-side sourcing is not available
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


# worker entry point that uses only the explicitly exported run context
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


# initialize a PSOCK cluster once without exporting the whole global environment
# Workers either source the exact same study scripts as the master session, or
# receive only the helper functions that are actually needed for the replication work.
make_simulation_cluster <- function(
    cores,
    script_dir = NULL,
    script_files = NULL,
    helper_env = environment(run_simulation_model_set)
) {

  cl <- parallel::makeCluster(cores)

  # cap implicit threading inside every worker before heavy work starts
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OMP_DYNAMIC = "FALSE",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      MKL_DYNAMIC = "FALSE",
      BLIS_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
    }

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

  # source the exact same worker script set when available
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

    # fallback for sessions where the scripts were sourced ad hoc and no project script_dir exists
    export_simulation_helpers_to_cluster(
      cl = cl,
      helper_env = helper_env
    )
  }

  cl
}


# export only the objects required by the current run instead of the entire global environment
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


# tune Elastic Net only once per unique Elastic Net tuning recipe and assign the result back to every model spec
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

  # if no Elastic Net model is present, nothing to do
  if (length(enet_specs) == 0) {
    return(specs)
  }

  # do we need to run internal tuning at all?
  needs_internal_tuning <- any(vapply(
    enet_specs,
    function(x) is.null(x$enet_tuning) && isTRUE(x$tune_enet),
    logical(1)
  ))

  df_pilot <- NULL

  if (needs_internal_tuning) {

    message("Stage 1.5/2: tuning Elastic Net")

    # simulate one pilot data set under the same data-generating mechanism
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

  # resolve one tuning object per unique tuning group
  tuning_group_ids <- sort(unique(vapply(enet_specs, function(x) x$enet_tuning_group_id, integer(1))))

  for (gid in tuning_group_ids) {

    idx <- which(vapply(specs, function(x) identical(x$enet_tuning_group_id, gid), logical(1)))
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

  # cap threaded math libraries in the master R session as well
  set_single_thread_math_env()
  cap_threads_runtime()

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

  # resolve Elastic Net tuning objects before we build the stage-1 groups
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

    # if available, use the exact project script set that the master session sourced
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

  # match choices
  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # create one model specification and hand it to the shared engine
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
