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


# determine the effective confounder order actually used in the full pipeline
compute_effective_c_order <- function(residualizer, sem_model, confounder_order) {

  # if no stage in the pipeline uses observed confounders, store 0
  if (residualizer == "none" && sem_model %in% c("riclpm", "dpm")) {
    return(0L)
  }

  as.integer(confounder_order)
}


# build an empty T-row result frame for a failed model within a replication
make_failed_replication_frame <- function(
    R,
    T,
    model_name,
    residualizer,
    sem_model,
    exclude,
    confounder_order,
    free_loadings,
    bootstrap_B,
    Phi,
    flag = 1L
) {

  # true lagged parameters are still known from Phi
  true_par <- extract_true_lagged_parameters(Phi = Phi, T = T)

  # convert the main-analysis flag into one-hot proportions
  flag_props <- flag_to_props(flag)

  data.frame(
    model_name = rep(model_name, T),
    R = rep(as.integer(R), T),
    T = seq_len(T),

    # analysis_flag retains the classification of the main fitted model itself
    analysis_flag = rep(as.integer(flag), T),

    # these three columns always sum to 1
    flag0 = rep(as.numeric(flag_props$flag0), T),
    flag1 = rep(as.numeric(flag_props$flag1), T),
    flag2 = rep(as.numeric(flag_props$flag2), T),

    # diagnostic column for the main fit
    improper_reason = rep(NA_character_, T),

    model = rep(encode_sem_model(sem_model), T),
    residualizer = rep(encode_residualizer(residualizer), T),
    exclusion = rep(collapse_exclusion(exclude), T),
    c_order = rep(compute_effective_c_order(residualizer, sem_model, confounder_order), T),
    free_loadings = rep(if (sem_model %in% c("riclpm", "dpm")) as.integer(free_loadings) else NA_integer_, T),
    bootstrap_B = rep(if (uses_bootstrap_se(residualizer)) as.integer(bootstrap_B) else 0L, T),
    bootstrap_prop_success = rep(NA_real_, T),
    BIC = rep(NA_real_, T),
    beta_x = true_par$beta_x,
    beta_y = true_par$beta_y,
    gamma_xy = true_par$gamma_xy,
    gamma_yx = true_par$gamma_yx,
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
    bootstrap_out = NULL
) {

  # classify the fit
  analysis_flag <- classify_fit_flag(fit_out$fit)

  # diagnose the main improper fit if applicable
  improper_reason <- if (analysis_flag == 2L) {
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

  # default bootstrap outputs
  bootstrap_prop_success <- NA_real_
  bootstrap_issue_vector <- NA_character_

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
      flag2 = bootstrap_out$flag2
    )
  }

  # true lagged parameters implied by Phi
  true_par <- extract_true_lagged_parameters(Phi = Phi, T = T)

  data.frame(
    model_name = rep(model_name, T),
    R = rep(as.integer(R), T),
    T = seq_len(T),

    # keep the main-analysis flag explicitly
    analysis_flag = rep(as.integer(analysis_flag), T),

    # these three columns always sum to 1
    flag0 = rep(as.numeric(flag_props$flag0), T),
    flag1 = rep(as.numeric(flag_props$flag1), T),
    flag2 = rep(as.numeric(flag_props$flag2), T),

    # one extra column for the main-fit improper reason
    improper_reason = rep(improper_reason, T),

    model = rep(encode_sem_model(spec$sem_model), T),
    residualizer = rep(encode_residualizer(spec$residualizer), T),
    exclusion = rep(collapse_exclusion(spec$exclude), T),
    c_order = rep(compute_effective_c_order(spec$residualizer, spec$sem_model, spec$confounder_order), T),
    free_loadings = rep(if (spec$sem_model %in% c("riclpm", "dpm")) as.integer(spec$free_loadings) else NA_integer_, T),
    bootstrap_B = rep(if (uses_bootstrap_se(spec$residualizer)) as.integer(spec$bootstrap_B) else 0L, T),
    bootstrap_prop_success = rep(bootstrap_prop_success, T),
    BIC = rep(extract_bic(fit_out$fit), T),
    beta_x = true_par$beta_x,
    beta_y = true_par$beta_y,
    gamma_xy = true_par$gamma_xy,
    gamma_yx = true_par$gamma_yx,
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
        exclude = spec$exclude,
        confounder_order = spec$confounder_order,
        free_loadings = spec$free_loadings,
        bootstrap_B = spec$bootstrap_B,
        Phi = Phi,
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
      confounder_order = proto$confounder_order,
      exclude = proto$exclude,
      xgb_tuning = proto$xgb_tuning,
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
        model_string = NULL
      )
    } else {
      fit_results[[spec$name]] <- fit_sem_on_prepared_data(
        df_prepared = prep$data,
        T = T,
        k = k,
        residualizer = spec$residualizer,
        sem_model = spec$sem_model,
        confounder_order = spec$confounder_order,
        exclude = spec$exclude,
        free_loadings = spec$free_loadings
      )
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
      bootstrap_out = bootstrap_out[[spec$name]]
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
    residualizer = c("none", "linear", "xgb"),
    sem_model = c("clpm", "riclpm", "dpm"),
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE,
    bootstrap_B = 50,
    bootstrap_seed = NULL,
    xgb_tuning = NULL,
    residualizer_args = list(),
    seed = NULL
) {

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  spec <- make_model_spec(
    name = "model_1",
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = bootstrap_B,
    xgb_tuning = xgb_tuning,
    tune_xgb = FALSE,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args
  )

  specs <- normalize_model_spec_list(list(spec))
  specs <- assign_xgb_tuning_group_ids(specs)
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