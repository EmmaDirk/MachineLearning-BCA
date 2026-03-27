# 09_one_replication_wrapper.R
# This wrapper executes one full replication of the simplified simulation study.
#
# One replication now means:
# 1) simulate one data set under one user-supplied Delta trajectory
# 2) fit exactly one analysis pipeline to that data set
# 3) if requested, replace model-based SEs by bootstrap SEs for two-stage methods
# 4) return a T-row data frame with one row per occasion
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


# build an empty T-row result frame for a failed replication
make_failed_replication_frame <- function(
    R,
    T,
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

  data.frame(
    R = rep(as.integer(R), T),
    T = seq_len(T),
    flag = rep(as.integer(flag), T),
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
    se_CYX = rep(NA_real_, T)
  )
}


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

  # match user choices
  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

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

  # if simulation failed, return a failed T-row frame immediately
  if (inherits(sim, "sim_error")) {
    return(make_failed_replication_frame(
      R = R,
      T = T,
      residualizer = residualizer,
      sem_model = sem_model,
      exclude = exclude,
      confounder_order = confounder_order,
      free_loadings = free_loadings,
      bootstrap_B = bootstrap_B,
      Phi = Phi,
      flag = 1L
    ))
  }

  # newer simulator may return a list with a data component; older versions may return the data directly
  df <- if (is.list(sim) && !is.data.frame(sim) && !is.null(sim$data)) sim$data else sim

  # fit the chosen analysis pipeline
  fit_out <- fit_analysis_pipeline(
    df = df,
    T = T,
    k = k,
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    xgb_tuning = xgb_tuning,
    residualizer_args = residualizer_args
  )

  # classify the fit
  flag <- classify_fit_flag(fit_out$fit)

  # extract model-based point estimates and standard errors
  lag <- extract_lagged_estimates(
    fit = fit_out$fit,
    T = T,
    model_type = sem_model
  )

  # default bootstrap success fraction
  bootstrap_prop_success <- NA_real_

  # replace standard errors by bootstrap SEs for two-stage methods if requested
  if (uses_bootstrap_se(residualizer) && bootstrap_B >= 2 && !is.null(fit_out$fit)) {

    se_boot <- bootstrap_pipeline_se(
      df = df,
      T = T,
      k = k,
      residualizer = residualizer,
      sem_model = sem_model,
      confounder_order = confounder_order,
      exclude = exclude,
      free_loadings = free_loadings,
      xgb_tuning = xgb_tuning,
      residualizer_args = residualizer_args,
      B = bootstrap_B,
      seed = bootstrap_seed
    )

    lag$se_ARX <- se_boot$ARX
    lag$se_ARY <- se_boot$ARY
    lag$se_CXY <- se_boot$CXY
    lag$se_CYX <- se_boot$CYX
    bootstrap_prop_success <- se_boot$bootstrap_prop_success
  }

  # true lagged parameters implied by Phi
  true_par <- extract_true_lagged_parameters(Phi = Phi, T = T)

  # build the requested final output for this replication
  data.frame(
    R = rep(as.integer(R), T),
    T = seq_len(T),
    flag = rep(as.integer(flag), T),
    model = rep(encode_sem_model(sem_model), T),
    residualizer = rep(encode_residualizer(residualizer), T),
    exclusion = rep(collapse_exclusion(exclude), T),
    c_order = rep(compute_effective_c_order(residualizer, sem_model, confounder_order), T),
    free_loadings = rep(if (sem_model %in% c("riclpm", "dpm")) as.integer(free_loadings) else NA_integer_, T),
    bootstrap_B = rep(if (uses_bootstrap_se(residualizer)) as.integer(bootstrap_B) else 0L, T),
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
    se_CYX = lag$se_CYX
  )
}