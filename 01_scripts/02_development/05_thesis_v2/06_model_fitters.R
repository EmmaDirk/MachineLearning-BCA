# These functions fit one or more analysis pipelines per simulated data set.
#
# The logic is now split into two clear stages:
# 1) prepare the stage-1 analysis data
#    - rename interaction columns for lavaan compatibility
#    - optionally residualise X and Y
# 2) fit the chosen SEM to that prepared data
#
# This split is important for efficiency:
# - multiple SEMs can now reuse the exact same prepared data
# - bootstrap refits can do the same within each resample
#
# The supported choices are:
# - residualizer = "none", "linear", or "xgb"
# - sem_model    = "clpm", "riclpm", or "dpm"
#
# The same analyst-side confounder choices are used throughout:
# - confounder_order says how rich the observed confounder set is
# - exclude removes observed confounders from that set
#
# Important:
# - CLPM can directly include observed confounders, but only when residualizer = "none"
# - RI-CLPM and DPM never directly include confounders here
# - when residualizer is "linear" or "xgb", the SEM is fit to residualised X and Y
# -------------------------------------------------------------------------------------------------

# rename interaction columns from c1:c2 to c1.c2 for lavaan model strings
rename_feature_columns_for_lavaan <- function(df) {

  # work on a data frame copy
  df <- as.data.frame(df)

  # lavaan model strings use dots in interaction names
  names(df) <- gsub(":", ".", names(df), fixed = TRUE)
  df
}


# rename excluded feature names the same way before they enter lavaan strings
rename_exclude_for_lavaan <- function(exclude) {

  if (is.null(exclude)) return(NULL)
  gsub(":", ".", exclude, fixed = TRUE)
}


# fit one lavaan model safely without stopping the full simulation
safe_fit_lavaan <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit the model
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

  # return fit and error
  list(fit = fit, err = err)
}


# validate the residualiser-specific argument set before calling stage 1
validate_residualizer_call <- function(
    residualizer,
    residualizer_args,
    xgb_tuning
) {

  # match the residualiser choice
  residualizer <- match.arg(residualizer, c("none", "linear", "xgb"))

  # normalize NULL to an empty list
  if (is.null(residualizer_args)) {
    residualizer_args <- list()
  }

  # require a list
  if (!is.list(residualizer_args)) {
    stop("'residualizer_args' must be a list.")
  }

  # allowed extra arguments by residualiser
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
    "nthread",
    "seed"
  )

  # no residualisation: no extra arguments should be supplied
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

  # linear residualisation: reject XGB-style arguments immediately
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

    return(invisible(TRUE))
  }

  # xgb residualisation: check allowed arguments and require tuning
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

  invisible(TRUE)
}


# apply the chosen residualiser once to the data set
apply_residualizer <- function(
    df,
    residualizer = c("none", "linear", "xgb"),
    k,
    confounder_order = 1,
    exclude = NULL,
    xgb_tuning = NULL,
    residualizer_args = list()
) {

  # match the residualiser choice
  residualizer <- match.arg(residualizer)

  # validate the residualiser-specific call before dispatch
  validate_residualizer_call(
    residualizer = residualizer,
    residualizer_args = residualizer_args,
    xgb_tuning = xgb_tuning
  )

  # no residualisation: return the data unchanged
  if (residualizer == "none") {
    return(df)
  }

  # linear residualisation
  if (residualizer == "linear") {

    # combine the shared analyst-side settings with any extra arguments
    args <- c(
      list(
        df = df,
        k = k,
        exclude = exclude,
        interaction_order = confounder_order
      ),
      residualizer_args
    )

    return(do.call(residualise_panel_linearC, args))
  }

  # XGB residualisation
  args <- c(
    list(
      df = df,
      tuning = xgb_tuning,
      k = k,
      exclude = exclude,
      interaction_order = confounder_order
    ),
    residualizer_args
  )

  do.call(residualise_panel_xgb, args)
}


# build the correct SEM model string from the chosen SEM option
build_sem_model_string <- function(
    T,
    sem_model = c("clpm", "riclpm", "dpm"),
    residualizer = c("none", "linear", "xgb"),
    k,
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE
) {

  # match arguments
  sem_model <- match.arg(sem_model)
  residualizer <- match.arg(residualizer)

  # CLPM can directly include observed confounders only when we did not residualise first
  if (sem_model == "clpm") {

    # after residualisation, fit a plain CLPM to the residualised series
    if (residualizer != "none") {
      return(build_clpm(T = T, k = 0, confounder_order = 0, exclude = NULL))
    }

    # otherwise fit the analyst-specified direct-adjustment CLPM
    return(build_clpm(
      T = T,
      k = k,
      confounder_order = confounder_order,
      exclude = rename_exclude_for_lavaan(exclude)
    ))
  }

  # RI-CLPM never directly includes confounders here
  if (sem_model == "riclpm") {
    return(build_riclpm(T = T, free_loadings = free_loadings))
  }

  # DPM never directly includes confounders here
  build_dpm(T = T, free_loadings = free_loadings)
}


# prepare the stage-1 analysis data once
prepare_analysis_data <- function(
    df,
    k,
    residualizer = c("none", "linear", "xgb"),
    confounder_order = 1,
    exclude = NULL,
    xgb_tuning = NULL,
    residualizer_args = list()
) {

  # match choices
  residualizer <- match.arg(residualizer)

  # harmonize data names once up front for any later direct-adjustment CLPM
  df_work <- rename_feature_columns_for_lavaan(df)

  # apply stage 1 safely
  df_stage2 <- tryCatch(
    apply_residualizer(
      df = df_work,
      residualizer = residualizer,
      k = k,
      confounder_order = confounder_order,
      exclude = rename_exclude_for_lavaan(exclude),
      xgb_tuning = xgb_tuning,
      residualizer_args = residualizer_args
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "stage1_error")
  )

  # if residualisation failed, return a structured failure result
  if (inherits(df_stage2, "stage1_error")) {
    return(list(
      data = NULL,
      err = df_stage2$message
    ))
  }

  list(
    data = df_stage2,
    err = NA_character_
  )
}


# fit the chosen SEM to already prepared stage-1 data
fit_sem_on_prepared_data <- function(
    df_prepared,
    T,
    k,
    residualizer = c("none", "linear", "xgb"),
    sem_model = c("clpm", "riclpm", "dpm"),
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE
) {

  # match choices
  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # if stage 1 failed upstream, return a structured failure fit result
  if (is.null(df_prepared)) {
    return(list(
      fit = NULL,
      err = "Prepared data is NULL.",
      data_used = NULL,
      model_string = NULL
    ))
  }

  # build the stage-2 model string
  model_string <- build_sem_model_string(
    T = T,
    sem_model = sem_model,
    residualizer = residualizer,
    k = k,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings
  )

  # fit the SEM safely
  fit_out <- safe_fit_lavaan(model_string = model_string, data = df_prepared)

  # return everything needed downstream
  list(
    fit = fit_out$fit,
    err = fit_out$err,
    data_used = df_prepared,
    model_string = model_string
  )
}


# fit one complete analysis pipeline to one data set
fit_analysis_pipeline <- function(
    df,
    T,
    k,
    residualizer = c("none", "linear", "xgb"),
    sem_model = c("clpm", "riclpm", "dpm"),
    confounder_order = 1,
    exclude = NULL,
    free_loadings = FALSE,
    xgb_tuning = NULL,
    residualizer_args = list()
) {

  # match choices
  residualizer <- match.arg(residualizer)
  sem_model <- match.arg(sem_model)

  # stage 1
  prep <- prepare_analysis_data(
    df = df,
    k = k,
    residualizer = residualizer,
    confounder_order = confounder_order,
    exclude = exclude,
    xgb_tuning = xgb_tuning,
    residualizer_args = residualizer_args
  )

  # if residualisation failed, stop here and return a structured failure result
  if (is.null(prep$data)) {
    return(list(
      fit = NULL,
      err = prep$err,
      data_used = NULL,
      model_string = NULL
    ))
  }

  # fit the SEM on the prepared data
  fit_sem_on_prepared_data(
    df_prepared = prep$data,
    T = T,
    k = k,
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings
  )
}
