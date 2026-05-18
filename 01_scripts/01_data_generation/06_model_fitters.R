# =================================================================================================
#
# These functions fit one or more analysis pipelines per simulated data set.
#
# The logic is split into two clear stages:
# 1) prepare the stage-1 analysis data
#    - rename interaction columns for lavaan compatibility
#    - optionally residualise X and Y
# 2) fit the chosen SEM to that prepared data
#
# The supported choices are:
# - residualizer = "none", "linear", "xgb", or "enet"
# - sem_model    = "clpm", "riclpm", or "dpm"
#
# The analyst-side confounder choices are now split into two separate layers:
# - sem_c_order / sem_exclude describe what the SEM layer directly uses
# - residualizer_c_order / residualizer_exclude describe what the residualiser sees
#
# Important:
# - CLPM, RI-CLPM, and DPM can directly include observed confounders, but only when residualizer = "none"
# - when residualizer is not "none", the SEM is fit to residualised X and Y
# =================================================================================================

# ---- lavaan-safe names --------------------------------------------------------------------------

# Rename interaction columns from Z1:2 to Z1.2 for lavaan model strings.

rename_feature_columns_for_lavaan <- function(df) {

  # Work on a data frame copy.

  df <- as.data.frame(df)

  # lavaan model strings use dots in interaction names.

  names(df) <- gsub(":", ".", names(df), fixed = TRUE)

  df
}


# Rename excluded feature names the same way before they enter lavaan strings.

rename_exclude_for_lavaan <- function(exclude) {

  if (is.null(exclude)) {
    return(NULL)
  }

  gsub(":", ".", exclude, fixed = TRUE)
}


# ---- lavaan fitting ------------------------------------------------------------------------------

# Fit one lavaan model safely without stopping the full simulation.

safe_fit_lavaan <- function(model_string, data) {

  # Initialize error message.

  err <- NA_character_

  # Try to fit the model.

  fit <- tryCatch(
    lavaan::lavaan(
      model     = model_string,
      data      = as.data.frame(data),
      estimator = "ML",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # Return fit and error.

  list(fit = fit, err = err)
}


# ---- panel structure helpers ---------------------------------------------------------------------

# Infer the number of observed waves from the prepared panel data.

infer_panel_T <- function(df, x_prefix = "x", y_prefix = "y") {

  # Work on a data frame copy.

  df <- as.data.frame(df)

  # Find all observed X and Y columns.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Extract their wave numbers.

  waves <- suppressWarnings(as.integer(sub("^[^0-9]+", "", c(x_cols, y_cols))))

  # If no waves were found, return 0.

  if (length(waves) == 0 || all(is.na(waves))) {
    return(0L)
  }

  as.integer(max(waves, na.rm = TRUE))
}


# ---- ML metric frames ----------------------------------------------------------------------------

# Build one empty T-row ML metric frame.

make_empty_ml_metric_frame <- function(T) {

  data.frame(
    T = seq_len(T),
    mse_x = rep(NA_real_, T),
    r2_x = rep(NA_real_, T),
    mse_y = rep(NA_real_, T),
    r2_y = rep(NA_real_, T),
    stringsAsFactors = FALSE
  )
}


# Standardize one ML metric frame to the requested set of T occasions.

standardize_ml_metric_frame <- function(ml_metrics, T) {

  # Start from the empty template.

  out <- make_empty_ml_metric_frame(T)

  # If no metrics were supplied, return the template.

  if (is.null(ml_metrics) || !is.data.frame(ml_metrics) || nrow(ml_metrics) == 0) {
    return(out)
  }

  # Require a time index.

  if (!("T" %in% names(ml_metrics))) {
    return(out)
  }

  # Keep only the supported columns.

  keep <- intersect(names(ml_metrics), names(out))
  tmp <- ml_metrics[, keep, drop = FALSE]

  # Merge onto the complete T grid.

  out <- merge(
    out["T"],
    tmp,
    by = "T",
    all.x = TRUE,
    sort = TRUE
  )

  # Ensure that all expected columns exist and are ordered correctly.

  for (nm in setdiff(names(make_empty_ml_metric_frame(T)), names(out))) {
    out[[nm]] <- NA_real_
  }

  out <- out[, names(make_empty_ml_metric_frame(T)), drop = FALSE]
  rownames(out) <- NULL

  out
}


# ---- residualiser validation ---------------------------------------------------------------------

# Validate the residualiser-specific argument set before calling stage 1.

validate_residualizer_call <- function(
    residualizer,
    residualizer_args,
    xgb_tuning,
    enet_tuning
) {

  # Match the residualiser choice.

  residualizer <- match.arg(residualizer, c("none", "linear", "xgb", "enet"))

  # Normalize NULL to an empty list.

  if (is.null(residualizer_args)) {
    residualizer_args <- list()
  }

  # Require a list.

  if (!is.list(residualizer_args)) {
    stop("'residualizer_args' must be a list.")
  }

  # Allowed extra arguments by residualiser.

  linear_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "seed"
  )

  xgb_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "use_oof",
    "late_start_wave",
    "nthread",
    "seed"
  )

  enet_allowed <- c(
    "x_prefix",
    "y_prefix",
    "c_prefix",
    "oof_folds",
    "seed"
  )

  # No residualisation: no extra arguments should be supplied.

  if (residualizer == "none") {

    if (length(residualizer_args) > 0) {
      stop(
        "You set residualizer = 'none', so 'residualizer_args' must be empty. ",
        "Unused entries: ",
        paste(names(residualizer_args), collapse = ", ")
      )
    }

    return(invisible(TRUE))
  }

  # Linear residualisation: reject XGB-style arguments immediately.

  if (residualizer == "linear") {

    bad <- setdiff(names(residualizer_args), linear_allowed)

    if (length(bad) > 0) {
      stop(
        "You set residualizer = 'linear', but these arguments are not valid for the linear residualiser: ",
        paste(bad, collapse = ", "),
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

    return(invisible(TRUE))
  }

  # XGB residualisation: check allowed arguments and require tuning.

  if (residualizer == "xgb") {

    bad <- setdiff(names(residualizer_args), xgb_allowed)

    if (length(bad) > 0) {
      stop(
        "You set residualizer = 'xgb', but these arguments are not valid for the XGB residualiser: ",
        paste(bad, collapse = ", "),
        ". Allowed extra arguments are: ",
        paste(xgb_allowed, collapse = ", "),
        "."
      )
    }

    if (is.null(xgb_tuning)) {
      stop("You set residualizer = 'xgb', so 'xgb_tuning' must be supplied.")
    }

    if (!is.null(enet_tuning)) {
      stop("You set residualizer = 'xgb', so 'enet_tuning' must be NULL.")
    }

    return(invisible(TRUE))
  }

  # Elastic Net residualisation: check allowed arguments and require tuning.

  if (residualizer == "enet") {

    bad <- setdiff(names(residualizer_args), enet_allowed)

    if (length(bad) > 0) {
      stop(
        "You set residualizer = 'enet', but these arguments are not valid for the Elastic Net residualiser: ",
        paste(bad, collapse = ", "),
        ". Allowed extra arguments are: ",
        paste(enet_allowed, collapse = ", "),
        "."
      )
    }

    if (is.null(enet_tuning)) {
      stop("You set residualizer = 'enet', so 'enet_tuning' must be supplied.")
    }

    if (!is.null(xgb_tuning)) {
      stop("You set residualizer = 'enet', so 'xgb_tuning' must be NULL.")
    }

    return(invisible(TRUE))
  }

  invisible(TRUE)
}


# ---- residualiser dispatch -----------------------------------------------------------------------

# Apply the chosen residualiser once to the data set.

apply_residualizer <- function(
    df,
    residualizer = c("none", "linear", "xgb", "enet"),
    k,
    residualizer_c_order = 1,
    residualizer_exclude = NULL,
    xgb_tuning = NULL,
    enet_tuning = NULL,
    residualizer_args = list()
) {

  # Match the residualiser choice.

  residualizer <- match.arg(residualizer)

  # Validate the residualiser-specific call before dispatch.

  validate_residualizer_call(
    residualizer = residualizer,
    residualizer_args = residualizer_args,
    xgb_tuning = xgb_tuning,
    enet_tuning = enet_tuning
  )

  # Allow custom prefixes when we need to build an empty metric frame.

  x_prefix <- if (!is.null(residualizer_args$x_prefix)) residualizer_args$x_prefix else "x"
  y_prefix <- if (!is.null(residualizer_args$y_prefix)) residualizer_args$y_prefix else "y"

  # Infer the number of observed waves once.

  T_panel <- infer_panel_T(df, x_prefix = x_prefix, y_prefix = y_prefix)

  # No residualisation: return the data unchanged together with empty ML metrics.

  if (residualizer == "none") {
    return(list(
      data = df,
      ml_metrics = make_empty_ml_metric_frame(T_panel)
    ))
  }

  # Linear residualisation.

  if (residualizer == "linear") {

    # Combine the shared analyst-side settings with any extra arguments.

    args <- c(
      list(
        df = df,
        k = k,
        exclude = residualizer_exclude,
        interaction_order = residualizer_c_order
      ),
      residualizer_args
    )

    return(do.call(residualise_panel_linearC, args))
  }

  # XGB residualisation.

  if (residualizer == "xgb") {

    args <- c(
      list(
        df = df,
        tuning = xgb_tuning,
        k = k,
        exclude = residualizer_exclude,
        interaction_order = residualizer_c_order
      ),
      residualizer_args
    )

    return(do.call(residualise_panel_xgb, args))
  }

  # Elastic Net residualisation.

  if (residualizer == "enet") {

    args <- c(
      list(
        df = df,
        tuning = enet_tuning,
        k = k,
        exclude = residualizer_exclude,
        interaction_order = residualizer_c_order
      ),
      residualizer_args
    )

    return(do.call(residualise_panel_enet, args))
  }
}


# ---- SEM model-string dispatch ------------------------------------------------------------------

# Build the correct SEM model string from the chosen SEM option.

build_sem_model_string <- function(
    T,
    sem_model = c("clpm", "riclpm", "dpm"),
    residualizer = c("none", "linear", "xgb", "enet"),
    k,
    sem_c_order = 1,
    sem_exclude = NULL,
    free_loadings = FALSE
) {

  # Match arguments.

  sem_model <- match.arg(sem_model)
  residualizer <- match.arg(residualizer)

  # CLPM can directly include observed confounders only when we did not residualise first.

  if (sem_model == "clpm") {

    # After residualisation, fit a plain CLPM to the residualised series.

    if (residualizer != "none") {
      return(build_clpm(T = T, k = 0, confounder_order = 0, exclude = NULL))
    }

    # Otherwise fit the analyst-specified direct-adjustment CLPM.

    return(build_clpm(
      T = T,
      k = k,
      confounder_order = sem_c_order,
      exclude = rename_exclude_for_lavaan(sem_exclude)
    ))
  }

  # After residualisation, fit the plain RI-CLPM to the residualised series.

  if (sem_model == "riclpm") {

    if (residualizer != "none") {
      return(build_riclpm(T = T, free_loadings = free_loadings, k = 0, confounder_order = 0, exclude = NULL))
    }

    # Otherwise fit the analyst-specified direct-adjustment RI-CLPM.

    return(build_riclpm(
      T = T,
      free_loadings = free_loadings,
      k = k,
      confounder_order = sem_c_order,
      exclude = rename_exclude_for_lavaan(sem_exclude)
    ))
  }

  # After residualisation, fit the plain DPM to the residualised series.

  if (residualizer != "none") {
    return(build_dpm(T = T, free_loadings = free_loadings, k = 0, confounder_order = 0, exclude = NULL))
  }

  # Otherwise fit the analyst-specified direct-adjustment DPM.

  build_dpm(
    T = T,
    free_loadings = free_loadings,
    k = k,
    confounder_order = sem_c_order,
    exclude = rename_exclude_for_lavaan(sem_exclude)
  )
}


# ---- stage-1 data preparation --------------------------------------------------------------------

# Prepare the stage-1 analysis data once.

prepare_analysis_data <- function(
    df,
    k,
    residualizer = c("none", "linear", "xgb", "enet"),
    residualizer_c_order = 1,
    residualizer_exclude = NULL,
    xgb_tuning = NULL,
    enet_tuning = NULL,
    residualizer_args = list()
) {

  # Match choices.

  residualizer <- match.arg(residualizer)

  # Allow custom prefixes when we need to build an empty metric frame.

  x_prefix <- if (!is.null(residualizer_args$x_prefix)) residualizer_args$x_prefix else "x"
  y_prefix <- if (!is.null(residualizer_args$y_prefix)) residualizer_args$y_prefix else "y"

  # Harmonize data names once up front for any later direct-adjustment CLPM.

  df_work <- rename_feature_columns_for_lavaan(df)

  # Infer the number of observed waves once.

  T_panel <- infer_panel_T(df_work, x_prefix = x_prefix, y_prefix = y_prefix)

  # Apply stage 1 safely.

  stage1_out <- tryCatch(
    apply_residualizer(
      df = df_work,
      residualizer = residualizer,
      k = k,
      residualizer_c_order = residualizer_c_order,
      residualizer_exclude = rename_exclude_for_lavaan(residualizer_exclude),
      xgb_tuning = xgb_tuning,
      enet_tuning = enet_tuning,
      residualizer_args = residualizer_args
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "stage1_error")
  )

  # If residualisation failed, return a structured failure result.

  if (inherits(stage1_out, "stage1_error")) {
    return(list(
      data = NULL,
      err = stage1_out$message,
      ml_metrics = make_empty_ml_metric_frame(T_panel)
    ))
  }

  list(
    data = stage1_out$data,
    err = NA_character_,
    ml_metrics = standardize_ml_metric_frame(stage1_out$ml_metrics, T = T_panel)
  )
}


# ---- stage-2 SEM fit ----------------------------------------------------------------------------

# Fit the chosen SEM to already prepared stage-1 data.

fit_sem_on_prepared_data <- function(
    df_prepared,
    T,
    k,
    residualizer = c("none", "linear", "xgb", "enet"),
    sem_model = c("clpm", "riclpm", "dpm"),
    sem_c_order = 1,
    sem_exclude = NULL,
    free_loadings = FALSE
) {

  # Match choices.

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # If stage 1 failed upstream, return a structured failure fit result.

  if (is.null(df_prepared)) {
    return(list(
      fit = NULL,
      err = "Prepared data is NULL.",
      data_used = NULL,
      model_string = NULL
    ))
  }

  # Build the stage-2 model string.

  model_string <- build_sem_model_string(
    T = T,
    sem_model = sem_model,
    residualizer = residualizer,
    k = k,
    sem_c_order = sem_c_order,
    sem_exclude = sem_exclude,
    free_loadings = free_loadings
  )

  # Fit the SEM safely.

  fit_out <- safe_fit_lavaan(model_string = model_string, data = df_prepared)

  # Return everything needed downstream.

  list(
    fit = fit_out$fit,
    err = fit_out$err,
    data_used = df_prepared,
    model_string = model_string
  )
}


# ---- full analysis pipeline ----------------------------------------------------------------------

# Fit one complete analysis pipeline to one data set.

fit_analysis_pipeline <- function(
    df,
    T,
    k,
    residualizer = c("none", "linear", "xgb", "enet"),
    sem_model = c("clpm", "riclpm", "dpm"),
    sem_c_order = 1,
    sem_exclude = NULL,
    residualizer_c_order = 1,
    residualizer_exclude = NULL,
    free_loadings = FALSE,
    xgb_tuning = NULL,
    enet_tuning = NULL,
    residualizer_args = list()
) {

  # Match choices.

  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # Stage 1.

  prep <- prepare_analysis_data(
    df = df,
    k = k,
    residualizer = residualizer,
    residualizer_c_order = residualizer_c_order,
    residualizer_exclude = residualizer_exclude,
    xgb_tuning = xgb_tuning,
    enet_tuning = enet_tuning,
    residualizer_args = residualizer_args
  )

  # If residualisation failed, stop here and return a structured failure result.

  if (is.null(prep$data)) {
    return(list(
      fit = NULL,
      err = prep$err,
      data_used = NULL,
      model_string = NULL,
      ml_metrics = prep$ml_metrics
    ))
  }

  # Fit the SEM on the prepared data.

  out <- fit_sem_on_prepared_data(
    df_prepared = prep$data,
    T = T,
    k = k,
    residualizer = residualizer,
    sem_model = sem_model,
    sem_c_order = sem_c_order,
    sem_exclude = sem_exclude,
    free_loadings = free_loadings
  )

  # Pass the stage-1 ML metrics through for optional downstream inspection.

  out$ml_metrics <- prep$ml_metrics

  out
}