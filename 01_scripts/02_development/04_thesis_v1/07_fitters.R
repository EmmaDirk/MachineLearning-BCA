# This script contains the direct fitting functions the user actually calls.
#
# The central design choice here is simplicity:
# - each function can be called directly on a data frame
# - each function returns only the small output that is usually needed
# - the common low-level work lives in 06_fit_helpers.R
#
# The three main user-facing functions are:
# - fit_clpm()
# - fit_riclpm()
# - fit_dpm()
#
# They all return a list with:
# - model_type
# - residualiser
# - converged
# - proper
# - bic
# - se_type
# - parameters
# - preprocess_info
# - error
#
# The returned parameter table is intentionally small:
# - one row per lagged parameter
# - estimate
# - standard error
# - optional model-based and bootstrap-based standard errors side by side
# ------------------------------------------------------------------------------------------

# CLPM
fit_clpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  # variable naming
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # which confounders the residualiser sees
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # which confounders the CLPM itself sees directly
  model_k = 0,
  model_exclude = NULL,
  model_confounder_order = 0,

  # bootstrap controls
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb controls
  xgb_tuning = NULL,
  xgb_tune_if_missing = TRUE,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning-stage controls
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan controls
  estimator = "ML",
  fixed_x = FALSE,

  # storage controls
  keep_syntax = FALSE,
  keep_lavaan_fit = FALSE
) {

  residualiser <- match.arg(residualiser)

  validate_common_fit_inputs(
    df = df,
    T = T,
    residualiser = residualiser,
    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,
    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,
    bootstrap_R = bootstrap_R,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds,
    estimator = estimator,
    fixed_x = fixed_x
  )

  validate_clpm_inputs(
    model_k = model_k,
    model_exclude = model_exclude,
    model_confounder_order = model_confounder_order
  )

  syntax <- build_model_syntax(
    model_type = "clpm",
    T = T,
    free_loadings = FALSE,
    model_k = model_k,
    model_exclude = model_exclude,
    model_confounder_order = model_confounder_order
  )

  fit_once_fun <- function(data_in) {
    fit_model_once(
      df = data_in,
      model_type = "clpm",
      T = T,
      syntax = syntax,
      residualiser = residualiser,
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix,
      resid_k = resid_k,
      resid_exclude = resid_exclude,
      resid_interaction_order = resid_interaction_order,
      xgb_tuning = xgb_tuning,
      xgb_tune_if_missing = xgb_tune_if_missing,
      xgb_tuning_grid = xgb_tuning_grid,
      xgb_cv_folds = xgb_cv_folds,
      xgb_nrounds_max = xgb_nrounds_max,
      xgb_early_stopping_rounds = xgb_early_stopping_rounds,
      xgb_oof_folds = xgb_oof_folds,
      xgb_nthread = xgb_nthread,
      xgb_seed = xgb_seed,
      estimator = estimator,
      fixed_x = fixed_x,
      keep_syntax = keep_syntax,
      keep_lavaan_fit = keep_lavaan_fit
    )
  }

  # if bootstrap is not requested, fit once and return
  if (is.null(bootstrap_R)) {
    return(fit_once_fun(df))
  }

  # otherwise use the bootstrap
  bootstrap_fit(
    df = df,
    fit_once_fun = fit_once_fun,
    R = bootstrap_R,
    seed = bootstrap_seed
  )
}

# RI-CLPM
fit_riclpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  # variable naming
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # which confounders the residualiser sees
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # measurement part
  free_loadings = FALSE,

  # bootstrap controls
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb controls
  xgb_tuning = NULL,
  xgb_tune_if_missing = TRUE,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning-stage controls
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan controls
  estimator = "ML",
  fixed_x = FALSE,

  # storage controls
  keep_syntax = FALSE,
  keep_lavaan_fit = FALSE
) {

  residualiser <- match.arg(residualiser)

  validate_common_fit_inputs(
    df = df,
    T = T,
    residualiser = residualiser,
    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,
    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,
    bootstrap_R = bootstrap_R,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds,
    estimator = estimator,
    fixed_x = fixed_x
  )

  validate_free_loading_input(free_loadings)

  syntax <- build_model_syntax(
    model_type = "riclpm",
    T = T,
    free_loadings = free_loadings
  )

  fit_once_fun <- function(data_in) {
    fit_model_once(
      df = data_in,
      model_type = "riclpm",
      T = T,
      syntax = syntax,
      residualiser = residualiser,
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix,
      resid_k = resid_k,
      resid_exclude = resid_exclude,
      resid_interaction_order = resid_interaction_order,
      xgb_tuning = xgb_tuning,
      xgb_tune_if_missing = xgb_tune_if_missing,
      xgb_tuning_grid = xgb_tuning_grid,
      xgb_cv_folds = xgb_cv_folds,
      xgb_nrounds_max = xgb_nrounds_max,
      xgb_early_stopping_rounds = xgb_early_stopping_rounds,
      xgb_oof_folds = xgb_oof_folds,
      xgb_nthread = xgb_nthread,
      xgb_seed = xgb_seed,
      estimator = estimator,
      fixed_x = fixed_x,
      keep_syntax = keep_syntax,
      keep_lavaan_fit = keep_lavaan_fit
    )
  }

  if (is.null(bootstrap_R)) {
    return(fit_once_fun(df))
  }

  bootstrap_fit(
    df = df,
    fit_once_fun = fit_once_fun,
    R = bootstrap_R,
    seed = bootstrap_seed
  )
}

# DPM
fit_dpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  # variable naming
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # which confounders the residualiser sees
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # measurement part
  free_loadings = FALSE,

  # bootstrap controls
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb controls
  xgb_tuning = NULL,
  xgb_tune_if_missing = TRUE,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning-stage controls
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan controls
  estimator = "ML",
  fixed_x = FALSE,

  # storage controls
  keep_syntax = FALSE,
  keep_lavaan_fit = FALSE
) {

  residualiser <- match.arg(residualiser)

  validate_common_fit_inputs(
    df = df,
    T = T,
    residualiser = residualiser,
    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,
    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,
    bootstrap_R = bootstrap_R,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds,
    estimator = estimator,
    fixed_x = fixed_x
  )

  validate_free_loading_input(free_loadings)

  syntax <- build_model_syntax(
    model_type = "dpm",
    T = T,
    free_loadings = free_loadings
  )

  fit_once_fun <- function(data_in) {
    fit_model_once(
      df = data_in,
      model_type = "dpm",
      T = T,
      syntax = syntax,
      residualiser = residualiser,
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix,
      resid_k = resid_k,
      resid_exclude = resid_exclude,
      resid_interaction_order = resid_interaction_order,
      xgb_tuning = xgb_tuning,
      xgb_tune_if_missing = xgb_tune_if_missing,
      xgb_tuning_grid = xgb_tuning_grid,
      xgb_cv_folds = xgb_cv_folds,
      xgb_nrounds_max = xgb_nrounds_max,
      xgb_early_stopping_rounds = xgb_early_stopping_rounds,
      xgb_oof_folds = xgb_oof_folds,
      xgb_nthread = xgb_nthread,
      xgb_seed = xgb_seed,
      estimator = estimator,
      fixed_x = fixed_x,
      keep_syntax = keep_syntax,
      keep_lavaan_fit = keep_lavaan_fit
    )
  }

  if (is.null(bootstrap_R)) {
    return(fit_once_fun(df))
  }

  bootstrap_fit(
    df = df,
    fit_once_fun = fit_once_fun,
    R = bootstrap_R,
    seed = bootstrap_seed
  )
}
