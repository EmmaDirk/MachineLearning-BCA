# 10_one_replication_wrapper.R
# This file executes one full replication for either:
# - one user-supplied model specification, or
# - a whole user-supplied set of model specifications
#
# One efficient replication means:
# 1) simulate one data set under one user-supplied Delta trajectory
# 2) prepare each unique stage-1 data set only once
# 3) fit every compatible SEM to that prepared data
# 4) if requested, bootstrap once per draw and reuse that draw across all models
# 5) return one stacked T-row data frame with one row per occasion per model
#
# The output data frame intentionally repeats study metadata across occasions.
# This makes later row-binding across model runs much simpler.
# -------------------------------------------------------------------------------------------------

# compact codes used in the final results frame
encode_sem_model <- function(sem_model) {

  switch(
    sem_model,
    clpm   = "C",
    riclpm = "R",
    dpm    = "D",
    stop("Unknown sem_model.")
  )
}


encode_residualizer <- function(residualizer) {

  switch(
    residualizer,
    none   = "N",
    linear = "L",
    xgb    = "X",
    enet   = "E",
    stop("Unknown residualizer.")
  )
}


# collapse the exclusion vector into one merge-friendly string
collapse_exclusion <- function(exclude) {

  if (is.null(exclude) || length(exclude) == 0) {
    return(NA_character_)
  }

  paste(exclude, collapse = "|")
}


# compute the metadata that are actually operative in the fitted pipeline
compute_effective_sem_exclusion <- function(residualizer, sem_model, sem_exclude) {

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))
  sem_model <- match.arg(sem_model, c("clpm", "riclpm", "dpm"))

  # after stage-1 residualisation, the SEM no longer receives observed confounders directly
  if (residualizer != "none") {
    return(NA_character_)
  }

  collapse_exclusion(sem_exclude)
}


compute_effective_sem_c_order <- function(residualizer, sem_model, sem_c_order) {

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))
  sem_model <- match.arg(sem_model, c("clpm", "riclpm", "dpm"))

  # after stage-1 residualisation, the SEM is fit without direct observed confounders
  if (residualizer != "none") {
    return(0L)
  }

  as.integer(sem_c_order[1])
}


compute_effective_residualizer_exclusion <- function(residualizer, residualizer_exclude) {

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))

  # without stage-1 residualisation, residualiser-side exclusions are not operative
  if (residualizer == "none") {
    return(NA_character_)
  }

  collapse_exclusion(residualizer_exclude)
}


compute_effective_residualizer_c_order <- function(residualizer, residualizer_c_order) {

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))

  # without stage-1 residualisation, no residualiser-side confounder columns are used
  if (residualizer == "none") {
    return(0L)
  }

  as.integer(residualizer_c_order[1])
}


# build an empty T-row result frame for a failed model within a replication
make_failed_replication_frame <- function(
    R,
    T,
    model_name,
    residualizer,
    sem_model,
    sem_exclude,
    sem_c_order,
    residualizer_exclude,
    residualizer_c_order,
    free_loadings,
    bootstrap_B,
    Phi,
    true_confounder_r2,
    flag = 1L
) {

  # true lagged parameters are still known from Phi
  true_par <- extract_true_lagged_parameters(Phi = Phi, T = T)

  # keep only the observed-wave true DGM quantities.
  # These are wave-specific and identical across fitted models within the same
  # replication because they are implied entirely by the DGM, not by the fitted model.
  true_r2_obs <- true_confounder_r2$observed_waves

  # convert the main-analysis flag into one-hot proportions
  flag_props <- flag_to_props(flag)

  data.frame(
    model_name = rep(as.character(model_name), T),
    R = rep(as.integer(R), T),
    T = seq_len(T),

    # analysis_flag retains the classification of the main fitted model itself
    analysis_flag = rep(as.integer(flag), T),

    # these four columns always sum to 1
    flag0 = rep(as.numeric(flag_props$flag0), T),
    flag1 = rep(as.numeric(flag_props$flag1), T),
    flag2 = rep(as.numeric(flag_props$flag2), T),
    flag3 = rep(as.numeric(flag_props$flag3), T),

    # diagnostic column for the main fit
    improper_reason = rep(NA_character_, T),

    model = rep(encode_sem_model(sem_model), T),
    residualizer = rep(encode_residualizer(residualizer), T),

    # explicit layer-specific metadata
    sem_exclusion = rep(compute_effective_sem_exclusion(residualizer, sem_model, sem_exclude), T),
    sem_c_order = rep(compute_effective_sem_c_order(residualizer, sem_model, sem_c_order), T),
    residualizer_exclusion = rep(compute_effective_residualizer_exclusion(residualizer, residualizer_exclude), T),
    residualizer_c_order = rep(compute_effective_residualizer_c_order(residualizer, residualizer_c_order), T),

    free_loadings = rep(if (sem_model %in% c("riclpm", "dpm")) as.integer(free_loadings) else NA_integer_, T),
    bootstrap_B = rep(if (uses_bootstrap_se(residualizer)) as.integer(bootstrap_B) else 0L, T),
    bootstrap_prop_success = rep(NA_real_, T),
    BIC = rep(NA_real_, T),
    beta_x = true_par$beta_x,
    beta_y = true_par$beta_y,
    gamma_xy = true_par$gamma_xy,
    gamma_yx = true_par$gamma_yx,
    internal_wave = true_r2_obs$t_internal,
    true_r2_x = true_r2_obs$true_r2_x,
    true_r2_y = true_r2_obs$true_r2_y,

    # store the direct true Delta_t coefficients as one list-column entry per wave.
    # Each element is a named numeric vector produced upstream by
    # compute_true_confounder_r2(). Because the DGM is shared across fitted models
    # within one replication, this column is repeated identically across models.
    true_delta_t_vector = I(true_r2_obs$true_delta_t_vector),
    mse_x = rep(NA_real_, T),
    r2_x = rep(NA_real_, T),
    mse_y = rep(NA_real_, T),
    r2_y = rep(NA_real_, T),
    ARX = rep(NA_real_, T),
    se_ARX = rep(NA_real_, T),
    ARY = rep(NA_real_, T),
    se_ARY = rep(NA_real_, T),
    CXY = rep(NA_real_, T),
    se_CXY = rep(NA_real_, T),
    CYX = rep(NA_real_, T),
    se_CYX = rep(NA_real_, T),
    bootstrap_issue_vector = I(rep(list(NA_character_), T)),
    stringsAsFactors = FALSE
  )
}


# collapse one fitted model into the final T-row output format
build_model_result_frame <- function(
    R,
    T,
    model_name,
    spec,
    Phi,
    fit_out,
    main_ml_metrics = NULL,
    bootstrap_out = NULL,
    true_confounder_r2
) {

  # classify the fit
  analysis_flag <- classify_fit_flag(fit_out$fit)

  # diagnose the main improper fit if applicable
  improper_reason <- if (analysis_flag %in% c(2L, 3L)) {
    diagnose_improper_fit(fit_out$fit)
  } else {
    NA_character_
  }

  # extract model-based point estimates and standard errors
  lag <- extract_lagged_estimates(
    fit = fit_out$fit,
    T = T,
    model_type = spec$sem_model
  )

  # standardize the main-run ML metrics to the final T-row output
  ml_main <- standardize_ml_metric_frame(main_ml_metrics, T = T)

  # default bootstrap outputs
  bootstrap_prop_success <- NA_real_
  bootstrap_issue_vector <- NA_character_

  # We do not store bootstrap summaries of the stage-1 OOF MSE / R^2 metrics
  # anymore. The final saved output focuses on the main per-replication OOF
  # point estimates, while Monte Carlo uncertainty is computed later across
  # replications. Bootstrap is still retained here for SEM-path standard errors
  # and bootstrap success-rate summaries.

  # default proportional flags are the one-hot encoding of the main fit
  flag_props <- flag_to_props(analysis_flag)

  # overwrite SEs and flag proportions when bootstrap output is available
  if (!is.null(bootstrap_out)) {

    lag$se_ARX <- bootstrap_out$ARX
    lag$se_ARY <- bootstrap_out$ARY
    lag$se_CXY <- bootstrap_out$CXY
    lag$se_CYX <- bootstrap_out$CYX

    bootstrap_prop_success <- bootstrap_out$bootstrap_prop_success
    bootstrap_issue_vector <- bootstrap_out$bootstrap_issue_vector

    flag_props <- list(
      flag0 = bootstrap_out$flag0,
      flag1 = bootstrap_out$flag1,
      flag2 = bootstrap_out$flag2,
      flag3 = bootstrap_out$flag3
    )

  }

  # true lagged parameters implied by Phi
  true_par <- extract_true_lagged_parameters(Phi = Phi, T = T)

  # keep only the observed-wave true DGM quantities.
  # This is the population benchmark implied by the DGM itself. It does not depend
  # on the fitted model, only on Phi, Delta_list, Omega11, Sigma, and burn-in.
  true_r2_obs <- true_confounder_r2$observed_waves

  data.frame(
    model_name = rep(model_name, T),
    R = rep(as.integer(R), T),
    T = seq_len(T),

    # keep the main-analysis flag explicitly
    analysis_flag = rep(as.integer(analysis_flag), T),

    # these four columns always sum to 1
    flag0 = rep(as.numeric(flag_props$flag0), T),
    flag1 = rep(as.numeric(flag_props$flag1), T),
    flag2 = rep(as.numeric(flag_props$flag2), T),
    flag3 = rep(as.numeric(flag_props$flag3), T),

    # one extra column for the main-fit improper reason
    improper_reason = rep(improper_reason, T),

    model = rep(encode_sem_model(spec$sem_model), T),
    residualizer = rep(encode_residualizer(spec$residualizer), T),

    # explicit layer-specific metadata
    sem_exclusion = rep(compute_effective_sem_exclusion(spec$residualizer, spec$sem_model, spec$sem_exclude), T),
    sem_c_order = rep(compute_effective_sem_c_order(spec$residualizer, spec$sem_model, spec$sem_c_order), T),
    residualizer_exclusion = rep(compute_effective_residualizer_exclusion(spec$residualizer, spec$residualizer_exclude), T),
    residualizer_c_order = rep(compute_effective_residualizer_c_order(spec$residualizer, spec$residualizer_c_order), T),

    free_loadings = rep(if (spec$sem_model %in% c("riclpm", "dpm")) as.integer(spec$free_loadings) else NA_integer_, T),
    bootstrap_B = rep(if (uses_bootstrap_se(spec$residualizer)) as.integer(spec$bootstrap_B) else 0L, T),
    bootstrap_prop_success = rep(bootstrap_prop_success, T),
    BIC = rep(extract_bic(fit_out$fit), T),
    beta_x = true_par$beta_x,
    beta_y = true_par$beta_y,
    gamma_xy = true_par$gamma_xy,
    gamma_yx = true_par$gamma_yx,
    internal_wave = true_r2_obs$t_internal,
    true_r2_x = true_r2_obs$true_r2_x,
    true_r2_y = true_r2_obs$true_r2_y,

    # keep the direct true Delta_t coefficients as well.
    #
    # This is intentionally a list-column rather than a wide set of scalar columns,
    # because the number of confounder features depends on the scenario and on
    # whether interactions are present. Storing one named vector per wave keeps the
    # final output compact while preserving the full information needed for later
    # inspection on the server output.
    true_delta_t_vector = I(true_r2_obs$true_delta_t_vector),
    mse_x = ml_main$mse_x,
    r2_x = ml_main$r2_x,
    mse_y = ml_main$mse_y,
    r2_y = ml_main$r2_y,
    ARX = lag$ARX,
    se_ARX = lag$se_ARX,
    ARY = lag$ARY,
    se_ARY = lag$se_ARY,
    CXY = lag$CXY,
    se_CXY = lag$se_CXY,
    CYX = lag$CYX,
    se_CYX = lag$se_CYX,
    bootstrap_issue_vector = I(rep(list(bootstrap_issue_vector), T)),
    stringsAsFactors = FALSE
  )
}


# run one efficient replication for a whole set of models
run_one_replication_model_set <- function(
    R,
    N,
    T,
    k,
    Phi,
    Sigma,
    Omega11,
    Delta_list,
    model_specs,
    stage1_groups,
    burn_in = 0L,
    bootstrap_seed = NULL,
    seed = NULL
) {

  # compute the true confounder trajectory object once for this replication
  # before any model is fit.
  #
  # This object is a property of the DGM itself and is therefore shared by all
  # fitted models within the same replication. It now contains both:
  # - the true confounder R^2 trajectory for X and Y
  # - the direct Delta_t coefficients for each wave in flattened vector form
  #
  # Computing it once here keeps the code efficient and ensures that these true
  # DGM quantities are still written to the output even if a later model fit fails.
  true_confounder_r2 <- compute_true_confounder_r2(
    T = T,
    Phi = Phi,
    Delta_list = Delta_list,
    Omega11 = Omega11,
    Sigma = Sigma,
    burn_in = burn_in
  )

  # simulate one data set safely
  sim <- tryCatch(
    simulate_panel_data(
      N = N,
      T = T,
      Phi = Phi,
      Delta_list = Delta_list,
      Omega11 = Omega11,
      Sigma = Sigma,
      burn_in = burn_in,
      seed = seed
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "sim_error")
  )

  # if simulation failed, return one failed T-row frame per model immediately
  if (inherits(sim, "sim_error")) {
    failed_list <- lapply(model_specs, function(spec) {
      make_failed_replication_frame(
        R = R,
        T = T,
        model_name = spec$name,
        residualizer = spec$residualizer,
        sem_model = spec$sem_model,
        sem_exclude = spec$sem_exclude,
        sem_c_order = spec$sem_c_order,
        residualizer_exclude = spec$residualizer_exclude,
        residualizer_c_order = spec$residualizer_c_order,
        free_loadings = spec$free_loadings,
        bootstrap_B = spec$bootstrap_B,
        Phi = Phi,
        true_confounder_r2 = true_confounder_r2,
        flag = 1L
      )
    })

    return(do.call(rbind, failed_list))
  }

  # the simulator returns the analysis data frame directly
  df <- sim

  # prepare each unique stage-1 data set exactly once
  prepared_by_group <- list()

  for (group_obj in stage1_groups) {

    proto <- group_obj$prototype

    prepared_by_group[[as.character(group_obj$stage1_group_id)]] <- prepare_analysis_data(
      df = df,
      k = k,
      residualizer = proto$residualizer,
      residualizer_c_order = proto$residualizer_c_order,
      residualizer_exclude = proto$residualizer_exclude,
      xgb_tuning = proto$xgb_tuning,
      enet_tuning = proto$enet_tuning,
      residualizer_args = proto$residualizer_args
    )
  }

  # fit every requested SEM on its already prepared stage-1 data
  fit_results <- list()

  for (spec in model_specs) {

    prep <- prepared_by_group[[as.character(spec$stage1_group_id)]]

    if (is.null(prep$data)) {
      fit_results[[spec$name]] <- list(
        fit = NULL,
        err = prep$err,
        data_used = NULL,
        model_string = NULL,
        ml_metrics = prep$ml_metrics
      )
    } else {
      fit_results[[spec$name]] <- fit_sem_on_prepared_data(
        df_prepared = prep$data,
        T = T,
        k = k,
        residualizer = spec$residualizer,
        sem_model = spec$sem_model,
        sem_c_order = spec$sem_c_order,
        sem_exclude = spec$sem_exclude,
        free_loadings = spec$free_loadings
      )

      # attach the stage-1 ML metrics used by this fitted pipeline
      fit_results[[spec$name]]$ml_metrics <- prep$ml_metrics
    }
  }

  # bootstrap only the models that need it and only when the main fit exists
  bootstrap_eligible <- Filter(function(spec) {
    uses_bootstrap_se(spec$residualizer) &&
      spec$bootstrap_B >= 2L &&
      !is.null(fit_results[[spec$name]]$fit)
  }, model_specs)

  bootstrap_out <- if (length(bootstrap_eligible) > 0) {
    bootstrap_model_set(
      df = df,
      T = T,
      k = k,
      model_specs = bootstrap_eligible,
      stage1_groups = stage1_groups,
      seed = bootstrap_seed
    )
  } else {
    setNames(vector("list", 0), character(0))
  }

  # build the final stacked results frame
  result_list <- lapply(model_specs, function(spec) {
    build_model_result_frame(
      R = R,
      T = T,
      model_name = spec$name,
      spec = spec,
      Phi = Phi,
      fit_out = fit_results[[spec$name]],
      main_ml_metrics = fit_results[[spec$name]]$ml_metrics,
      bootstrap_out = bootstrap_out[[spec$name]],
      true_confounder_r2 = true_confounder_r2
    )
  })

  do.call(rbind, result_list)
}




# convenience wrapper for running one single model through the shared replication engine
run_one_replication <- function(
    R,
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
    xgb_tuning = NULL,
    enet_tuning = NULL,
    residualizer_args = list(),
    seed = NULL,
    tune_xgb = FALSE,
    tune_enet = FALSE,
    xgb_tune_args = list(),
    enet_tune_args = list()
) {

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

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

  specs <- normalize_model_spec_list(list(spec))
  validate_model_spec_list(specs)

  base_seed_for_tuning <- if (is.null(seed) || length(seed) == 0 || is.na(seed[1])) {
    1234L
  } else {
    as.integer(seed[1])
  }

  specs <- resolve_xgb_tuning_objects(
    model_specs = specs,
    N = N,
    T = T,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    burn_in = burn_in,
    base_seed = base_seed_for_tuning
  )

  specs <- resolve_enet_tuning_objects(
    model_specs = specs,
    N = N,
    T = T,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    burn_in = burn_in,
    base_seed = base_seed_for_tuning
  )

  specs <- assign_stage1_group_ids(specs)
  stage1_groups <- build_stage1_groups(specs)

  out <- run_one_replication_model_set(
    R = R,
    N = N,
    T = T,
    k = k,
    Phi = Phi,
    Sigma = Sigma,
    Omega11 = Omega11,
    Delta_list = Delta_list,
    model_specs = specs,
    stage1_groups = stage1_groups,
    burn_in = burn_in,
    bootstrap_seed = bootstrap_seed,
    seed = seed
  )

  out
}
