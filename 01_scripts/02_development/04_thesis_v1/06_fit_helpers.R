# This script contains the small helpers used by the direct fitting functions.
#
# The aim here is to keep the top-level fitters readable:
# - fit_clpm()
# - fit_riclpm()
# - fit_dpm()
#
# We therefore place the repeated low-level work here:
# 1) argument checks
# 2) confounder bookkeeping
# 3) model syntax selection
# 4) residualisation preparation
# 5) lavaan fitting
# 6) convergence / properness checks
# 7) parameter extraction
# 8) bootstrap standard errors for two-stage procedures
#
# The output is intentionally small.
# By default we only keep:
# - the lagged parameter table
# - a convergence flag
# - a proper-solution flag
# - the BIC
# - a small amount of metadata
# ------------------------------------------------------------------------------------------

# helper: use x unless it is NULL
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# helper: check whether something is a single whole number
is_whole_number_scalar <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x == as.integer(x)
}

# helper: safe wrapper around lavInspect()
# purpose:
# we often want to inspect a fitted lavaan object, but this can fail for bad fits.
# rather than stopping everything, we return NULL and let the diagnostics handle it.
safe_lavinspect <- function(fit, what) {
  tryCatch(
    lavaan::lavInspect(fit, what),
    error = function(e) NULL
  )
}

# helper: safely subset a parameter table, creating missing columns as NA
# purpose:
# lavaan sometimes omits columns when a fit is odd or incomplete.
# this helper keeps the extractor stable.
safe_parameter_estimates_subset <- function(pe, needed_cols) {
  for (nm in needed_cols) {
    if (!nm %in% names(pe)) {
      pe[[nm]] <- NA
    }
  }

  pe[, needed_cols, drop = FALSE]
}

# helper: validate the parts of the arguments that all fitters share
validate_common_fit_inputs <- function(
  df,
  T,
  residualiser,
  x_prefix,
  y_prefix,
  c_prefix,
  resid_k,
  resid_exclude,
  resid_interaction_order,
  bootstrap_R,
  xgb_oof_folds,
  xgb_nthread,
  xgb_cv_folds,
  xgb_nrounds_max,
  xgb_early_stopping_rounds,
  estimator,
  fixed_x
) {

  # data must be coercible to a data frame
  df <- as.data.frame(df)

  # T must be at least 2, because otherwise there are no lagged effects
  if (!is_whole_number_scalar(T) || T < 2) {
    stop("'T' must be a single integer >= 2.")
  }

  # residualiser must be one of the implemented options
  if (!(residualiser %in% c("none", "linear", "xgb"))) {
    stop("'residualiser' must be one of 'none', 'linear', or 'xgb'.")
  }

  # prefixes must be small non-empty strings
  for (val in list(x_prefix = x_prefix, y_prefix = y_prefix, c_prefix = c_prefix)) {
    if (!is.character(val[[1]]) || length(val[[1]]) != 1 || nchar(val[[1]]) == 0) {
      stop("Variable prefixes must be single non-empty character strings.")
    }
  }

  # interaction order controls which confounder terms the residualiser sees
  if (!(resid_interaction_order %in% c(1, 2, 3))) {
    stop("'resid_interaction_order' must be 1, 2, or 3.")
  }

  # if resid_k is supplied, it must be a positive integer
  if (!is.null(resid_k)) {
    if (!is_whole_number_scalar(resid_k) || resid_k < 1) {
      stop("'resid_k' must be NULL or a single integer >= 1.")
    }
  }

  # exclude vectors must be character vectors if supplied
  if (!is.null(resid_exclude) && !is.character(resid_exclude)) {
    stop("'resid_exclude' must be NULL or a character vector.")
  }

  # bootstrap count can be NULL, otherwise it must be positive
  if (!is.null(bootstrap_R)) {
    if (!is_whole_number_scalar(bootstrap_R) || bootstrap_R < 1) {
      stop("'bootstrap_R' must be NULL or a single integer >= 1.")
    }
  }

  # xgb cost controls must be valid, even if xgb is not used
  if (!is_whole_number_scalar(xgb_oof_folds) || xgb_oof_folds < 2) {
    stop("'xgb_oof_folds' must be a single integer >= 2.")
  }

  if (!is_whole_number_scalar(xgb_nthread) || xgb_nthread < 1) {
    stop("'xgb_nthread' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(xgb_cv_folds) || xgb_cv_folds < 2) {
    stop("'xgb_cv_folds' must be a single integer >= 2.")
  }

  if (!is_whole_number_scalar(xgb_nrounds_max) || xgb_nrounds_max < 1) {
    stop("'xgb_nrounds_max' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(xgb_early_stopping_rounds) ||
      xgb_early_stopping_rounds < 1) {
    stop("'xgb_early_stopping_rounds' must be a single integer >= 1.")
  }

  # lavaan controls
  if (!is.character(estimator) || length(estimator) != 1) {
    stop("'estimator' must be a single character string.")
  }

  if (!is.logical(fixed_x) || length(fixed_x) != 1 || is.na(fixed_x)) {
    stop("'fixed_x' must be TRUE or FALSE.")
  }

  # check that the data has the required x and y columns
  x_cols <- paste0(x_prefix, 1:T)
  y_cols <- paste0(y_prefix, 1:T)

  missing_xy <- setdiff(c(x_cols, y_cols), names(df))
  if (length(missing_xy) > 0) {
    stop(
      "The following required panel columns are missing: ",
      paste(missing_xy, collapse = ", ")
    )
  }

  invisible(df)
}

# helper: validate the extra CLPM arguments
validate_clpm_inputs <- function(model_k, model_exclude, model_confounder_order) {

  if (!is_whole_number_scalar(model_k) || model_k < 0) {
    stop("'model_k' must be a single integer >= 0.")
  }

  if (!is.null(model_exclude) && !is.character(model_exclude)) {
    stop("'model_exclude' must be NULL or a character vector.")
  }

  if (!(model_confounder_order %in% c(0, 1, 2, 3))) {
    stop("'model_confounder_order' must be 0, 1, 2, or 3.")
  }
}

# helper: validate the free-loading toggle for RI-CLPM and DPM
validate_free_loading_input <- function(free_loadings) {
  if (!is.logical(free_loadings) || length(free_loadings) != 1 || is.na(free_loadings)) {
    stop("'free_loadings' must be TRUE or FALSE.")
  }
}

# helper: choose the columns that the residualiser is allowed to see
# purpose:
# the user explicitly said they like careful control over which confounders are visible.
# this helper implements that control in one place.
get_used_confounders <- function(df, k = NULL, c_prefix = "c", exclude = NULL) {

  df <- as.data.frame(df)

  # if k is NULL, we take every base confounder column that matches the prefix
  if (is.null(k)) {
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)
  } else {
    c_cols <- paste0(c_prefix, 1:k)
  }

  if (length(c_cols) == 0) {
    stop("No confounder columns found.")
  }

  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0) {
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))
  }

  if (!is.null(exclude)) {
    missing_exclude <- setdiff(exclude, names(df))
    if (length(missing_exclude) > 0) {
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))
    }

    c_cols <- setdiff(c_cols, exclude)
  }

  if (length(c_cols) == 0) {
    stop("No confounders left after exclusion.")
  }

  c_cols
}

# helper: build a small metadata list describing the preprocessing step
build_preprocess_info <- function(method, confounders_used, exclude, interaction_order, tuned_now = FALSE) {
  list(
    residualiser = method,
    confounders_used = confounders_used,
    exclude = exclude,
    interaction_order = interaction_order,
    tuned_now = tuned_now
  )
}

# helper: rename interaction columns so lavaan can see them safely
# purpose:
# the simulated data use names like c1:c2 and c1:c2:c3.
# for model.matrix() this is convenient, but inside lavaan syntax colons are awkward.
# we therefore translate every colon to a dot right before fitting.
# example:
# - c1:c2   becomes c1.c2
# - c1:c2:c3 becomes c1.c2.c3
rename_columns_for_lavaan <- function(df) {

  df <- as.data.frame(df)
  names(df) <- gsub(":", ".", names(df), fixed = TRUE)
  df
}

# helper: residualise or leave the data untouched
# purpose:
# this keeps the model-fitting functions small.
prepare_analysis_data <- function(
  df,
  residualiser = "none",
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,
  xgb_tuning = NULL,
  xgb_tune_if_missing = TRUE,
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123
) {

  df <- as.data.frame(df)

  # no residualisation: simply return the data as is
  if (residualiser == "none") {
    return(list(
      data = df,
      preprocess_info = build_preprocess_info(
        method = "none",
        confounders_used = NULL,
        exclude = NULL,
        interaction_order = NULL,
        tuned_now = FALSE
      ),
      xgb_tuning = NULL
    ))
  }

  # identify which confounders the residualiser may use
  used_c <- get_used_confounders(
    df = df,
    k = resid_k,
    c_prefix = c_prefix,
    exclude = resid_exclude
  )

  # linear residualisation is the simple benchmark
  if (residualiser == "linear") {
    data_out <- residualise_panel_linearC(
      df = df,
      k = resid_k,
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix,
      exclude = resid_exclude,
      interaction_order = resid_interaction_order
    )

    return(list(
      data = data_out,
      preprocess_info = build_preprocess_info(
        method = "linear",
        confounders_used = used_c,
        exclude = resid_exclude,
        interaction_order = resid_interaction_order,
        tuned_now = FALSE
      ),
      xgb_tuning = NULL
    ))
  }

  # xgb residualisation is more flexible but also more expensive
  if (residualiser == "xgb") {

    tuned_now <- FALSE

    # if no tuning object was supplied, we can tune on the current data
    if (is.null(xgb_tuning)) {
      if (!isTRUE(xgb_tune_if_missing)) {
        stop("XGB residualisation requested, but no tuning object was supplied.")
      }

      xgb_tuning <- tune_residualise_panel_xgb(
        df = df,
        k = resid_k,
        x_prefix = x_prefix,
        y_prefix = y_prefix,
        c_prefix = c_prefix,
        exclude = resid_exclude,
        interaction_order = resid_interaction_order,
        tuning_grid = xgb_tuning_grid,
        cv_folds = xgb_cv_folds,
        nrounds_max = xgb_nrounds_max,
        early_stopping_rounds = xgb_early_stopping_rounds,
        nthread = xgb_nthread,
        seed = xgb_seed
      )

      tuned_now <- TRUE
    }

    data_out <- residualise_panel_xgb(
      df = df,
      tuning = xgb_tuning,
      k = resid_k,
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix,
      exclude = resid_exclude,
      interaction_order = resid_interaction_order,
      oof_folds = xgb_oof_folds,
      nthread = xgb_nthread,
      seed = xgb_seed
    )

    return(list(
      data = data_out,
      preprocess_info = build_preprocess_info(
        method = "xgb",
        confounders_used = used_c,
        exclude = resid_exclude,
        interaction_order = resid_interaction_order,
        tuned_now = tuned_now
      ),
      xgb_tuning = xgb_tuning
    ))
  }

  stop("Unknown residualiser.")
}

# helper: decide which model string to build
build_model_syntax <- function(
  model_type,
  T,
  free_loadings = FALSE,
  model_k = 0,
  model_exclude = NULL,
  model_confounder_order = 0
) {

  if (model_type == "clpm") {
    return(build_clpm(
      T = T,
      k = model_k,
      confounder_order = model_confounder_order,
      exclude = model_exclude
    ))
  }

  if (model_type == "riclpm") {
    return(build_riclpm(
      T = T,
      free_loadings = free_loadings
    ))
  }

  if (model_type == "dpm") {
    return(build_dpm(
      T = T,
      free_loadings = free_loadings
    ))
  }

  stop("Unsupported model type.")
}

# helper: run lavaan
run_lavaan_fit <- function(data, syntax, estimator = "ML", fixed_x = FALSE) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required.")
  }

  lavaan::sem(
    model = syntax,
    data = as.data.frame(data),
    estimator = estimator,
    fixed.x = fixed_x,
    meanstructure = FALSE
  )
}

# helper: detect whether the solution is improper
# purpose:
# convergence alone is not enough.
# a model may converge numerically and still be unusable.
detect_improper_solution <- function(fit) {

  converged <- tryCatch(
    lavaan::lavInspect(fit, "converged"),
    error = function(e) FALSE
  )

  # if there is no convergence, we do not treat the fit as proper
  if (!isTRUE(converged)) {
    return(TRUE)
  }

  pe <- tryCatch(
    lavaan::parameterEstimates(fit),
    error = function(e) NULL
  )

  if (is.null(pe)) {
    return(TRUE)
  }

  # check for negative variances
  variances <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]
  negative_variance <- any(variances$est < 0, na.rm = TRUE)

  # check whether the latent covariance matrix is non-positive definite
  cov_lv <- safe_lavinspect(fit, "cov.lv")
  non_pd_lv <- FALSE

  if (!is.null(cov_lv) && is.matrix(cov_lv)) {
    eig <- tryCatch(
      eigen(cov_lv, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    non_pd_lv <- any(is.na(eig)) || any(eig <= 0)
  }

  # check whether the residual covariance matrix is non-positive definite
  theta <- safe_lavinspect(fit, "theta")
  non_pd_theta <- FALSE

  if (!is.null(theta) && is.matrix(theta)) {
    eig <- tryCatch(
      eigen(theta, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    non_pd_theta <- any(is.na(eig)) || any(eig <= 0)
  }

  isTRUE(negative_variance || non_pd_lv || non_pd_theta)
}

# helper: extract the small diagnostics object we care about
extract_fit_diagnostics <- function(lavaan_fit) {

  converged <- tryCatch(
    lavaan::lavInspect(lavaan_fit, "converged"),
    error = function(e) FALSE
  )

  bic <- tryCatch(
    stats::BIC(lavaan_fit),
    error = function(e) NA_real_
  )

  proper <- tryCatch(
    !detect_improper_solution(lavaan_fit),
    error = function(e) FALSE
  )

  list(
    converged = isTRUE(converged),
    proper = isTRUE(proper),
    bic = bic
  )
}

# helper: build the lagged paths that should be extracted
# purpose:
# this lets us return the same parameter table structure even when a fit fails.
build_target_paths <- function(model_type, T, x_prefix = "x", y_prefix = "y") {

  rows <- list()

  # CLPM and DPM express the lagged effects with the observed x / y variables
  if (model_type %in% c("clpm", "dpm")) {
    for (t in 2:T) {

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0(x_prefix, t),
        rhs = paste0(x_prefix, t - 1),
        path_type = "ar_x",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0(x_prefix, t),
        rhs = paste0(y_prefix, t - 1),
        path_type = "cl_y_to_x",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0(y_prefix, t),
        rhs = paste0(x_prefix, t - 1),
        path_type = "cl_x_to_y",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0(y_prefix, t),
        rhs = paste0(y_prefix, t - 1),
        path_type = "ar_y",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )
    }

    return(do.call(rbind, rows))
  }

  # RI-CLPM expresses the lagged effects with within-person latent variables
  if (model_type == "riclpm") {
    for (t in 2:T) {

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0("w", x_prefix, t),
        rhs = paste0("w", x_prefix, t - 1),
        path_type = "ar_x",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0("w", x_prefix, t),
        rhs = paste0("w", y_prefix, t - 1),
        path_type = "cl_y_to_x",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0("w", y_prefix, t),
        rhs = paste0("w", x_prefix, t - 1),
        path_type = "cl_x_to_y",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )

      rows[[length(rows) + 1]] <- data.frame(
        lhs = paste0("w", y_prefix, t),
        rhs = paste0("w", y_prefix, t - 1),
        path_type = "ar_y",
        wave_from = t - 1,
        wave_to = t,
        stringsAsFactors = FALSE
      )
    }

    return(do.call(rbind, rows))
  }

  stop("Unsupported model type.")
}

# helper: create an empty parameter table with the right shape
# purpose:
# when fitting fails, downstream code still gets the same table layout.
empty_parameter_table <- function(model_type, T, x_prefix = "x", y_prefix = "y") {

  out <- build_target_paths(
    model_type = model_type,
    T = T,
    x_prefix = x_prefix,
    y_prefix = y_prefix
  )

  out$param_id <- paste0(out$lhs, "~", out$rhs)
  out$est <- NA_real_
  out$se <- NA_real_
  out$se_model <- NA_real_
  out$se_boot <- NA_real_

  out[, c(
    "param_id", "lhs", "rhs", "path_type", "wave_from", "wave_to",
    "est", "se", "se_model", "se_boot"
  )]
}

# helper: extract only the lagged parameters we care about
extract_model_parameters <- function(lavaan_fit, model_type, T, x_prefix = "x", y_prefix = "y") {

  pe <- lavaan::parameterEstimates(lavaan_fit, standardized = FALSE)

  target <- build_target_paths(
    model_type = model_type,
    T = T,
    x_prefix = x_prefix,
    y_prefix = y_prefix
  )

  needed_cols <- c("lhs", "op", "rhs", "est", "se")
  pe_sub <- safe_parameter_estimates_subset(pe, needed_cols)

  out <- merge(
    target,
    pe_sub,
    by.x = c("lhs", "rhs"),
    by.y = c("lhs", "rhs"),
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[out$op == "~" | is.na(out$op), , drop = FALSE]
  out$param_id <- paste0(out$lhs, "~", out$rhs)
  out$se_model <- out$se
  out$se_boot <- NA_real_

  out <- out[, c(
    "param_id", "lhs", "rhs", "path_type", "wave_from", "wave_to",
    "est", "se", "se_model", "se_boot"
  )]

  rownames(out) <- NULL
  out
}

# helper: fit a model once and return the minimal output structure
fit_model_once <- function(
  df,
  model_type,
  T,
  syntax,
  residualiser = "none",
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,
  xgb_tuning = NULL,
  xgb_tune_if_missing = TRUE,
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,
  estimator = "ML",
  fixed_x = FALSE,
  keep_syntax = FALSE,
  keep_lavaan_fit = FALSE
) {

  # prepare the data first, but capture errors here as well
  # purpose:
  # fitting can fail before lavaan is even called, for example during
  # residualisation or xgb tuning. In that case we still want a small
  # failure object with the usual parameter-table shape.
  prepared <- tryCatch(
    prepare_analysis_data(
      df = df,
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
      xgb_seed = xgb_seed
    ),
    error = function(e) e
  )

  if (inherits(prepared, "error")) {
    out <- list(
      model_type = model_type,
      residualiser = residualiser,
      converged = FALSE,
      proper = FALSE,
      bic = NA_real_,
      se_type = "model",
      parameters = empty_parameter_table(
        model_type = model_type,
        T = T,
        x_prefix = x_prefix,
        y_prefix = y_prefix
      ),
      preprocess_info = NULL,
      error = conditionMessage(prepared)
    )

    if (isTRUE(keep_syntax)) {
      out$syntax <- syntax
    }

    if (isTRUE(keep_lavaan_fit)) {
      out$lavaan_fit <- NULL
    }

    return(out)
  }

  # before fitting, make the data-column names lavaan-safe
  data_for_fit <- rename_columns_for_lavaan(prepared$data)

  # fit the model, but capture errors
  fit <- tryCatch(
    run_lavaan_fit(
      data = data_for_fit,
      syntax = syntax,
      estimator = estimator,
      fixed_x = fixed_x
    ),
    error = function(e) e
  )

  # if fitting failed, return a small failure object
  if (inherits(fit, "error")) {
    out <- list(
      model_type = model_type,
      residualiser = residualiser,
      converged = FALSE,
      proper = FALSE,
      bic = NA_real_,
      se_type = "model",
      parameters = empty_parameter_table(
        model_type = model_type,
        T = T,
        x_prefix = x_prefix,
        y_prefix = y_prefix
      ),
      preprocess_info = prepared$preprocess_info,
      error = conditionMessage(fit)
    )

    if (isTRUE(keep_syntax)) {
      out$syntax <- syntax
    }

    if (isTRUE(keep_lavaan_fit)) {
      out$lavaan_fit <- NULL
    }

    if (residualiser == "xgb") {
      out$xgb_tuning <- prepared$xgb_tuning
    }

    return(out)
  }

  # extract diagnostics and parameter table
  diagnostics <- extract_fit_diagnostics(fit)
  params <- extract_model_parameters(
    lavaan_fit = fit,
    model_type = model_type,
    T = T,
    x_prefix = x_prefix,
    y_prefix = y_prefix
  )

  out <- list(
    model_type = model_type,
    residualiser = residualiser,
    converged = diagnostics$converged,
    proper = diagnostics$proper,
    bic = diagnostics$bic,
    se_type = "model",
    parameters = params,
    preprocess_info = prepared$preprocess_info,
    error = NA_character_
  )

  if (isTRUE(keep_syntax)) {
    out$syntax <- syntax
  }

  if (isTRUE(keep_lavaan_fit)) {
    out$lavaan_fit <- fit
  }

  if (residualiser == "xgb") {
    out$xgb_tuning <- prepared$xgb_tuning
  }

  out
}

# helper: draw one bootstrap sample
draw_bootstrap_sample <- function(df) {
  idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
  as.data.frame(df[idx, , drop = FALSE])
}

# helper: keep only the estimate column from a fitted result
extract_bootstrap_draw <- function(fit_result) {
  out <- fit_result$parameters[, c("param_id", "est"), drop = FALSE]
  rownames(out) <- NULL
  out
}

# helper: replace the model SEs by bootstrap SEs
summarise_bootstrap_parameters <- function(original_parameters, bootstrap_draws) {

  all_draws <- do.call(rbind, lapply(seq_along(bootstrap_draws), function(i) {

    d <- bootstrap_draws[[i]]

    if (is.null(d)) {
      return(NULL)
    }

    d$replicate <- i
    d
  }))

  # if every bootstrap fit failed, keep the structure but set se_boot to NA
  if (is.null(all_draws) || nrow(all_draws) == 0) {
    out <- original_parameters
    out$se_model <- out$se
    out$se_boot <- NA_real_
    out$se <- NA_real_
    return(out)
  }

  boot_se <- aggregate(est ~ param_id, data = all_draws, FUN = stats::sd, na.rm = TRUE)
  names(boot_se)[names(boot_se) == "est"] <- "se_boot"

  out <- merge(
    original_parameters,
    boot_se,
    by = "param_id",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_new")
  )

  # merge() can create a duplicate se_boot column name when the input already had one
  if ("se_boot_new" %in% names(out)) {
    out$se_boot <- out$se_boot_new
    out$se_boot_new <- NULL
  }

  out$se_model <- out$se
  out$se <- out$se_boot

  out <- out[, c(
    "param_id", "lhs", "rhs", "path_type", "wave_from", "wave_to",
    "est", "se", "se_model", "se_boot"
  )]

  rownames(out) <- NULL
  out
}

# helper: run the bootstrap for a fitted model function
# purpose:
# model-based SEs are fine for one-stage fits.
# for two-stage fits, the bootstrap lets us carry the uncertainty of both stages.
bootstrap_fit <- function(df, fit_once_fun, R, seed = 123) {

  set.seed(seed)

  # first fit on the original data
  original_fit <- fit_once_fun(df)

  # if the original fit itself failed, return it unchanged
  if (!isTRUE(original_fit$converged)) {
    original_fit$se_type <- "bootstrap"
    original_fit$parameters$se_model <- original_fit$parameters$se
    original_fit$parameters$se_boot <- NA_real_
    original_fit$parameters$se <- NA_real_
    original_fit$bootstrap_summary <- list(
      requested_replicates = R,
      successful_replicates = 0L
    )
    return(original_fit)
  }

  boot_draws <- vector("list", R)
  success <- logical(R)

  for (b in seq_len(R)) {

    boot_df <- draw_bootstrap_sample(df)

    fit_b <- tryCatch(
      fit_once_fun(boot_df),
      error = function(e) e
    )

    if (inherits(fit_b, "error")) {
      boot_draws[[b]] <- NULL
      success[b] <- FALSE
      next
    }

    # we keep the draw even when the fit is improper,
    # because the user said convergence and properness should be distinct.
    # for the bootstrap SD itself, the main issue is whether an estimate exists.
    boot_draws[[b]] <- extract_bootstrap_draw(fit_b)
    success[b] <- TRUE
  }

  original_fit$parameters <- summarise_bootstrap_parameters(
    original_parameters = original_fit$parameters,
    bootstrap_draws = boot_draws
  )

  original_fit$se_type <- "bootstrap"
  original_fit$bootstrap_summary <- list(
    requested_replicates = R,
    successful_replicates = sum(success)
  )

  original_fit
}
