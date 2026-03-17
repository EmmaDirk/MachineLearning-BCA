# the purpose of this script is to write a function that fits the lavaan models, with
# the option to residualise using linear models or XGBoost models.
# As such we get the following functions:
#
# 1) CLPM fitter with:
#     a) no residualisation
#     b) linear residualisation
#     c) XGBoost residualisation
#
# 2) DPM fitter with:
#     a) no residualisation
#     b) no residualisation, but free loadings
#     c) linear residualisation
#     d) XGBoost residualisation
#
# 3) RI-CLPM fitter with:
#     a) no residualisation
#     b) no residualisation, but free loadings
#     c) linear residualisation
#     d) XGBoost residualisation
#
# The functions devised here sometimes combine two models. In that case, uncertainty for both models
# needs to be taken into account. For this reason, we need bootstrap estimates for the s.e. parameters.
# To do this, we will here design a pipeline of functions:
#
# 1) Small helpers.
# 2) Function that determines which model to fit, and checks.
# 3) XGB tuner wrapper.
# 4) lavaan model wrapper.
# 5) residualiser wrapper.
# 6) lavaan fitter
# 7) fit checks
# 8) extractor for fit stats and parameters
# 9) fit wrapper
# 10) bootstrap setup + wrapper
# 11) top level wrapper
# 12) the functions to run each model

# -------------------------------- 1) Small helpers ----------------------------
# function that says: use X if not null, otherwise use Y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# check if x is a whole number
is_whole_number_scalar <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x == as.integer(x)
}

# use trycatch to avoid errors
safe_lavinspect <- function(fit, what) {
  tryCatch(
    lavaan::lavInspect(fit, what),
    error = function(e) NULL
  )
}

# helper that returns the parameter-estimate columns we need, creating missing ones as NA
safe_parameter_estimates_subset <- function(pe, needed_cols) {
  for (nm in needed_cols) {
    if (!nm %in% names(pe)) {
      pe[[nm]] <- NA
    }
  }

  pe[, needed_cols, drop = FALSE]
}

# helper to construct the target path table used by the extractor
build_target_paths <- function(spec) {

  model_type <- spec$model$model_type
  T <- spec$model$T
  x_prefix <- spec$variables$x_prefix
  y_prefix <- spec$variables$y_prefix

  rows <- list()

  # CLPM and DPM use observed x / y names
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

  # RI-CLPM uses within-person latent names wx / wy
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

# -------------------------------- 2) Argument Checks ----------------------------
make_sem_spec <- function(
  model_type,                                # "clpm", "riclpm", "dpm"
  T,                                         # number of waves
  residualiser = c("none", "linear", "xgb"), # "none", "linear", "xgb"

  # variable naming
  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # residualiser confounders
  resid_k = NULL,                            # number of confounders
  resid_exclude = NULL,                      # confounders the residualiser does not see
  resid_interaction_order = 1,               # order of interaction terms the residualiser sees

  # CLPM outcome-model confounders
  model_k = 0,                               # number of confounders
  model_exclude = NULL,                      # confounders the outcome-model does not see
  model_confounder_order = 0,                # order of interaction terms the outcome-model sees

  # feeing up CLPM / DPM loadings
  free_loadings = FALSE,

  # lavaan fitting
  estimator = "ML",
  fixed_x = FALSE,                           # allow exogenous variables to have variance and covariance

  # bootstrap
  bootstrap_R = NULL,                        # number of bootstrap samples
  bootstrap_seed = 123,

  # XGB fitting stage
  xgb_tuning = NULL,                         # pass the tuned model
  xgb_oof_folds = 2,                         # number of out-of-fold splits
  xgb_nthread = 1,                           # number of threads (set to 1 always in parallel simulations)
  xgb_seed = 123,

  # XGB tuning-stage controls
  xgb_tuning_grid = NULL,                    # pass the grid
  xgb_cv_folds = 5,                          # number of CV folds
  xgb_nrounds_max = 400,                     # maximum number of boosting iterations
  xgb_early_stopping_rounds = 20             # early stopping rounds
) {

  # argument checks
  residualiser <- match.arg(residualiser)

  analysis_stages <- if (residualiser == "none") 1L else 2L

  # check if bootstrap is required
  requires_bootstrap <- analysis_stages == 2L

  # save the passed arguments
  spec <- list(
    model = list(
      model_type = model_type,
      T = T,
      free_loadings = free_loadings
    ),

    variables = list(
      x_prefix = x_prefix,
      y_prefix = y_prefix,
      c_prefix = c_prefix
    ),

    residualisation = list(
      method = residualiser,
      k = resid_k,
      exclude = resid_exclude,
      interaction_order = resid_interaction_order
    ),

    outcome_model = list(
      k = model_k,
      exclude = model_exclude,
      confounder_order = model_confounder_order
    ),

    fitting = list(
      estimator = estimator,
      fixed_x = fixed_x
    ),

    bootstrap = list(
      requires_bootstrap = requires_bootstrap,
      R = bootstrap_R,
      seed = bootstrap_seed
    ),

    xgb = list(
      tuning = xgb_tuning,
      oof_folds = xgb_oof_folds,
      nthread = xgb_nthread,
      seed = xgb_seed,

      tuning_grid = xgb_tuning_grid,
      cv_folds = xgb_cv_folds,
      nrounds_max = xgb_nrounds_max,
      early_stopping_rounds = xgb_early_stopping_rounds
    ),

    pipeline = list(
      analysis_stages = analysis_stages
    )
  )

  spec
}

# function that checks the arguments
validate_sem_spec <- function(spec) {

  # ---------- model type ----------
  if (!(spec$model$model_type %in% c("clpm", "riclpm", "dpm"))) {
    stop("'model_type' must be one of 'clpm', 'riclpm', 'dpm'.")
  }

  # ---------- T ----------
  if (!is_whole_number_scalar(spec$model$T) || spec$model$T < 2) {
    stop("'T' must be a single integer >= 2.")
  }

  # ---------- prefixes ----------
  for (nm in c("x_prefix", "y_prefix", "c_prefix")) {
    val <- spec$variables[[nm]]
    if (!is.character(val) || length(val) != 1 || nchar(val) == 0) {
      stop(sprintf("'%s' must be a single non-empty character string.", nm))
    }
  }

  # ---------- residualiser ----------
  if (!(spec$residualisation$method %in% c("none", "linear", "xgb"))) {
    stop("'residualiser' must be one of 'none', 'linear', 'xgb'.")
  }

  # ---------- interaction orders ----------
  if (!(spec$residualisation$interaction_order %in% c(1, 2, 3))) {
    stop("'resid_interaction_order' must be 1, 2, or 3.")
  }

  if (!(spec$outcome_model$confounder_order %in% c(0, 1, 2, 3))) {
    stop("'model_confounder_order' must be 0, 1, 2, or 3.")
  }

  # ---------- k values ----------
  if (!is.null(spec$residualisation$k)) {
    if (!is_whole_number_scalar(spec$residualisation$k) || spec$residualisation$k < 1) {
      stop("'resid_k' must be NULL or a single integer >= 1.")
    }
  }

  if (!is_whole_number_scalar(spec$outcome_model$k) || spec$outcome_model$k < 0) {
    stop("'model_k' must be a single integer >= 0.")
  }

  # ---------- free loadings ----------
  if (isTRUE(spec$model$free_loadings) &&
      !(spec$model$model_type %in% c("riclpm", "dpm"))) {
    stop("'free_loadings = TRUE' is only allowed for 'riclpm' and 'dpm'.")
  }

  # ---------- CLPM direct confounders only ----------
  if (spec$model$model_type != "clpm") {
    if (!identical(spec$outcome_model$k, 0) ||
        !is.null(spec$outcome_model$exclude) ||
        !identical(spec$outcome_model$confounder_order, 0)) {
      warning("Direct outcome-model confounder settings are ignored for non-CLPM models.")
    }
  }

  # ---------- bootstrap required for two-stage ----------
  if (isTRUE(spec$bootstrap$requires_bootstrap)) {
    if (!is_whole_number_scalar(spec$bootstrap$R) || spec$bootstrap$R < 1) {
      stop("Two-stage analyses require bootstrap. Please provide 'bootstrap_R >= 1'.")
    }
  }

  # ---------- XGB fitting controls ----------
  if (!is_whole_number_scalar(spec$xgb$oof_folds) || spec$xgb$oof_folds < 2) {
    stop("'xgb_oof_folds' must be a single integer >= 2.")
  }

  if (!is_whole_number_scalar(spec$xgb$nthread) || spec$xgb$nthread < 1) {
    stop("'xgb_nthread' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(spec$xgb$cv_folds) || spec$xgb$cv_folds < 2) {
    stop("'xgb_cv_folds' must be a single integer >= 2.")
  }

  if (!is_whole_number_scalar(spec$xgb$nrounds_max) || spec$xgb$nrounds_max < 1) {
    stop("'xgb_nrounds_max' must be a single integer >= 1.")
  }

  if (!is_whole_number_scalar(spec$xgb$early_stopping_rounds) ||
      spec$xgb$early_stopping_rounds < 1) {
    stop("'xgb_early_stopping_rounds' must be a single integer >= 1.")
  }

  # ---------- xgb tuning presence warning ----------
  if (spec$residualisation$method == "xgb" && is.null(spec$xgb$tuning)) {
    warning(
      "XGB residualisation requested, but no tuning object is attached to spec yet. ",
      "This is fine before tuning. Before fitting, attach fixed tuning."
    )
  }

  spec
}

# -------------------------------- 3) XGB tuner wrapper -------------------------------
tune_xgb_once <- function(df, spec) {

  # check the arguments first
  spec <- validate_sem_spec(spec)

  # check that the residualiser is XGB
  if (spec$residualisation$method != "xgb") {
    stop("tune_xgb_once() should only be used when residualiser = 'xgb'.")
  }

  # run the XGB tuner
  tuning <- tune_residualise_panel_xgb(
    df = df,
    k = spec$residualisation$k,
    x_prefix = spec$variables$x_prefix,
    y_prefix = spec$variables$y_prefix,
    c_prefix = spec$variables$c_prefix,
    exclude = spec$residualisation$exclude,
    interaction_order = spec$residualisation$interaction_order,
    tuning_grid = spec$xgb$tuning_grid,
    cv_folds = spec$xgb$cv_folds,
    nrounds_max = spec$xgb$nrounds_max,
    early_stopping_rounds = spec$xgb$early_stopping_rounds,
    nthread = spec$xgb$nthread,
    seed = spec$xgb$seed
  )

  # attach metadata for compatibility checks
  tuning$.pipeline_meta <- list(
    x_prefix = spec$variables$x_prefix,
    y_prefix = spec$variables$y_prefix,
    c_prefix = spec$variables$c_prefix,
    resid_k = spec$residualisation$k,
    resid_exclude = spec$residualisation$exclude,
    resid_interaction_order = spec$residualisation$interaction_order,
    xgb_tuning_grid = spec$xgb$tuning_grid,
    xgb_cv_folds = spec$xgb$cv_folds,
    xgb_nrounds_max = spec$xgb$nrounds_max,
    xgb_early_stopping_rounds = spec$xgb$early_stopping_rounds,
    xgb_nthread = spec$xgb$nthread,
    xgb_seed = spec$xgb$seed
  )

  tuning
}

# check that the output of the tuning function is compatible
check_xgb_tuning_compatibility <- function(spec, tuning) {

  # check that a tuning object was supplied
  if (is.null(tuning)) {
    stop("No XGB tuning object supplied.")
  }

  # check that the tuning object is compatible
  meta <- tuning$.pipeline_meta
  if (is.null(meta)) {
    warning("Tuning object has no '.pipeline_meta'. Skipping strict compatibility check.")
    return(TRUE)
  }

  ok <- TRUE

  ok <- ok && identical(meta$x_prefix, spec$variables$x_prefix)
  ok <- ok && identical(meta$y_prefix, spec$variables$y_prefix)
  ok <- ok && identical(meta$c_prefix, spec$variables$c_prefix)
  ok <- ok && identical(meta$resid_k, spec$residualisation$k)
  ok <- ok && identical(sort(meta$resid_exclude %||% character(0)),
                        sort(spec$residualisation$exclude %||% character(0)))
  ok <- ok && identical(meta$resid_interaction_order, spec$residualisation$interaction_order)

  # if fails one of the checks, stop
  if (!ok) {
    stop("The supplied XGB tuning object is not compatible with the current residualiser specification.")
  }

  TRUE
}

# -------------------------------- 4) build lavaan model -------------------------------
build_model_syntax_from_spec <- function(spec) {

  # check the arguments
  spec <- validate_sem_spec(spec)

  # get the model type and T
  model_type <- spec$model$model_type
  T <- spec$model$T

  # if model type is clpm, use the clpm builder
  if (model_type == "clpm") {
    return(
      build_clpm(
        T = T,
        k = spec$outcome_model$k,
        confounder_order = spec$outcome_model$confounder_order,
        exclude = spec$outcome_model$exclude
      )
    )
  }

  # if model type is riclpm, use the riclpm builder
  if (model_type == "riclpm") {
    return(
      build_riclpm(
        T = T,
        free_loadings = spec$model$free_loadings
      )
    )
  }

  # if model type is dpm, use the dpm builder
  if (model_type == "dpm") {
    return(
      build_dpm(
        T = T,
        free_loadings = spec$model$free_loadings
      )
    )
  }

  stop("Unsupported model type.")
}

# -------------------------------- 5) residualisers --------------------------------------
# a function to get which confounders the residualiser sees.
get_used_confounders <- function(
  df,                                           # data frame
  k = NULL,                                     # number of confounders
  c_prefix = "c",                               # confounder prefix
  exclude = NULL) {                             # confounders to exclude

  # if k is null, use all confounders
  if (is.null(k)) {
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)
  } else {
    c_cols <- paste0(c_prefix, 1:k)
  }

  # stop if no confounders found
  if (length(c_cols) == 0) {
    stop("No confounder columns found.")
  }

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0) {
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))
  }

  # exclude confounders
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      stop("'exclude' must be a character vector.")
    }

    missing_exclude <- setdiff(exclude, names(df))
    if (length(missing_exclude) > 0) {
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))
    }

    c_cols <- setdiff(c_cols, exclude)
  }

  # stop if no confounders left
  if (length(c_cols) == 0) {
    stop("No confounders left after exclusion.")
  }

  c_cols
}

# wrapper around the linear residualiser using the specifications
residualise_data_linear <- function(df, spec) {

  out <- residualise_panel_linearC(
    df = df,
    k = spec$residualisation$k,
    x_prefix = spec$variables$x_prefix,
    y_prefix = spec$variables$y_prefix,
    c_prefix = spec$variables$c_prefix,
    exclude = spec$residualisation$exclude,
    interaction_order = spec$residualisation$interaction_order
  )

  used_c <- get_used_confounders(
    df = df,
    k = spec$residualisation$k,
    c_prefix = spec$variables$c_prefix,
    exclude = spec$residualisation$exclude
  )

  list(
    data = out,
    preprocess_info = list(
      residualiser = "linear",
      confounders_used = used_c,
      exclude = spec$residualisation$exclude,
      interaction_order = spec$residualisation$interaction_order
    )
  )
}

# wrapper around the XGB residualiser using the specifications
residualise_data_xgb <- function(df, spec) {

  if (is.null(spec$xgb$tuning)) {
    stop("XGB residualisation requested, but no fixed tuning object was supplied.")
  }

  check_xgb_tuning_compatibility(spec, spec$xgb$tuning)

  out <- residualise_panel_xgb(
    df = df,
    tuning = spec$xgb$tuning,
    k = spec$residualisation$k,
    x_prefix = spec$variables$x_prefix,
    y_prefix = spec$variables$y_prefix,
    c_prefix = spec$variables$c_prefix,
    exclude = spec$residualisation$exclude,
    interaction_order = spec$residualisation$interaction_order,
    oof_folds = spec$xgb$oof_folds,
    nthread = spec$xgb$nthread,
    seed = spec$xgb$seed
  )

  used_c <- get_used_confounders(
    df = df,
    k = spec$residualisation$k,
    c_prefix = spec$variables$c_prefix,
    exclude = spec$residualisation$exclude
  )

  list(
    data = out,
    preprocess_info = list(
      residualiser = "xgb",
      confounders_used = used_c,
      exclude = spec$residualisation$exclude,
      interaction_order = spec$residualisation$interaction_order,
      oof_folds = spec$xgb$oof_folds,
      tuning_grid = spec$xgb$tuning_grid,
      cv_folds = spec$xgb$cv_folds,
      nrounds_max = spec$xgb$nrounds_max,
      early_stopping_rounds = spec$xgb$early_stopping_rounds,
      tuning_summary = spec$xgb$tuning$final
    )
  )
}

# wrapper that calls the different residualiser functions
prepare_analysis_data <- function(df, spec) {

  method <- spec$residualisation$method

  if (method == "none") {
    return(
      list(
        data = as.data.frame(df),
        preprocess_info = list(
          residualiser = "none",
          confounders_used = NULL,
          exclude = NULL,
          interaction_order = NULL
        )
      )
    )
  }

  if (method == "linear") {
    return(residualise_data_linear(df, spec))
  }

  if (method == "xgb") {
    return(residualise_data_xgb(df, spec))
  }

  stop("Unknown residualiser.")
}

# ------------------------------ 6) lavaan fitter ------------------------------------------
run_lavaan_fit <- function(data, syntax, spec) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required.")
  }

  lavaan::sem(
    model = syntax,
    data = data,
    estimator = spec$fitting$estimator,
    fixed.x = spec$fitting$fixed_x,
    meanstructure = FALSE
  )
}

# ------------------------------ 7) fit checks -----------------------------------------
# detect if the fit is not properly identified
detect_improper_solution <- function(fit) {

  # check if the model converged
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) return(TRUE)

  # get the parameter estimates
  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)
  if (is.null(pe)) return(TRUE)

  # check for negative variances
  variances <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]
  negative_variance <- any(variances$est < 0, na.rm = TRUE)

  # check for non-positive definite latent variances
  cov_lv <- safe_lavinspect(fit, "cov.lv")
  non_pd_lv <- FALSE

  if (!is.null(cov_lv) && is.matrix(cov_lv)) {
    eig <- tryCatch(
      eigen(cov_lv, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    non_pd_lv <- any(is.na(eig)) || any(eig <= 0)
  }

  # check for non-positive definite residual covariance matrix
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

# ------------------------------- 8) extractors -----------------------------------------
# first get the BIC
extract_fit_diagnostics <- function(lavaan_fit) {

  # check if the model converged
  converged <- tryCatch(
    lavaan::lavInspect(lavaan_fit, "converged"),
    error = function(e) FALSE
  )

  # get the BIC
  bic <- tryCatch(stats::BIC(lavaan_fit), error = function(e) NA_real_)

  # detect if the model is not properly identified
  improper <- tryCatch(detect_improper_solution(lavaan_fit), error = function(e) TRUE)

  data.frame(
    converged = isTRUE(converged),
    improper = isTRUE(improper),
    bic = bic,
    stringsAsFactors = FALSE
  )
}

# now get the parameters of the autoregressive and cross-lagged paths
extract_sem_paths <- function(lavaan_fit, spec) {

  # get the parameter estimates
  pe <- lavaan::parameterEstimates(lavaan_fit, standardized = FALSE)

  # build the target path table for the requested model
  target <- build_target_paths(spec)

  # keep the columns we need; if lavaan omitted one (for instance "label"), create it as NA
  needed_cols <- c("lhs", "op", "rhs", "label", "est", "se", "z", "pvalue")
  pe_sub <- safe_parameter_estimates_subset(pe, needed_cols)

  # merge target paths with parameter estimates
  out <- merge(
    target,
    pe_sub,
    by.x = c("lhs", "rhs"),
    by.y = c("lhs", "rhs"),
    all.x = TRUE,
    sort = FALSE
  )

  # keep matched regression rows, but also retain unmatched target rows as NA
  out <- out[out$op == "~" | is.na(out$op), , drop = FALSE]
  out$param_id <- paste0(out$lhs, "~", out$rhs)

  out <- out[, c(
    "param_id", "lhs", "rhs", "path_type", "wave_from", "wave_to",
    "label", "est", "se", "z", "pvalue"
  )]

  rownames(out) <- NULL
  out
}

# now a top level wrapper bringing together the fit and path extractor.
extract_sem_results <- function(lavaan_fit, spec) {
  list(
    paths = extract_sem_paths(lavaan_fit, spec),
    diagnostics = extract_fit_diagnostics(lavaan_fit)
  )
}

# --------------------------------- 9) fit wrapper ---------------------------------------------
# top level wrapper to fit and extract results
fit_sem_once <- function(df, spec) {

  # first check the arguments
  spec <- validate_sem_spec(spec)

  # if residualisation method is XGB, check that the tuning is compatible
  if (spec$residualisation$method == "xgb") {
    check_xgb_tuning_compatibility(spec, spec$xgb$tuning)
  }

  # now call all functions
  # first residualise
  prepared <- prepare_analysis_data(df, spec)

  # get the lavaan syntax
  syntax <- build_model_syntax_from_spec(spec)

  # fit the model
  fit <- run_lavaan_fit(prepared$data, syntax, spec)

  # extract results
  extracted <- extract_sem_results(fit, spec)

  # put that all in a list
  list(
    spec = spec,
    syntax = syntax,
    preprocess_info = prepared$preprocess_info,
    data_used = prepared$data,
    lavaan_fit = fit,
    diagnostics = extracted$diagnostics,
    paths = extracted$paths
  )
}

# ------------------------------- 10) bootstrap setup ------------------------------------------
# function to draw a bootstrap sample
draw_bootstrap_sample <- function(df) {

  # draw a bootstrap sample as large as the original, but with replacement
  idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
  as.data.frame(df[idx, , drop = FALSE])
}

# function to extract bootstrap targets
extract_bootstrap_targets <- function(fit_object) {
  out <- fit_object$paths[, c("param_id", "est"), drop = FALSE]
  rownames(out) <- NULL
  out
}

# function to replace s.e. estimate with bootstrap estimate
summarise_bootstrap_paths <- function(original_paths, bootstrap_draws) {

  # combine all bootstrap draws
  all_draws <- do.call(rbind, lapply(seq_along(bootstrap_draws), function(i) {

    # get the bootstrap draw
    d <- bootstrap_draws[[i]]

    # if the draw is NULL, return NULL
    if (is.null(d)) return(NULL)

    # add a replicate column
    d$replicate <- i
    d
  }))

  # if something went wrong, store NA
  if (is.null(all_draws) || nrow(all_draws) == 0) {
    out <- original_paths
    out$se_model <- out$se
    out$se_boot <- NA_real_
    out$se <- NA_real_
    return(out)
  }

  # the bootstrap standard error
  boot_se <- aggregate(est ~ param_id, data = all_draws, FUN = stats::sd, na.rm = TRUE)
  names(boot_se)[names(boot_se) == "est"] <- "se_boot"

  # combine
  out <- merge(original_paths, boot_se, by = "param_id", all.x = TRUE, sort = FALSE)
  out$se_model <- out$se
  out$se <- out$se_boot

  # reorder
  out <- out[, c(
    "param_id", "lhs", "rhs", "path_type", "wave_from", "wave_to",
    "label", "est", "se_model", "se_boot", "se", "z", "pvalue"
  )]

  rownames(out) <- NULL
  out
}

# function to summarise bootstrap failures
summarise_bootstrap_failures <- function(rep_results, R) {

  # save the status: success, nonconvergence, improper
  status <- data.frame(
    replicate = seq_len(R),
    success = sapply(rep_results, function(x) isTRUE(x$success)),
    converged = sapply(rep_results, function(x) isTRUE(x$converged)),
    improper = sapply(rep_results, function(x) isTRUE(x$improper)),
    stringsAsFactors = FALSE
  )

  list(
    requested_replicates = R,
    successful_replicates = sum(status$success),
    failed_replicates = sum(!status$success),
    nonconverged_replicates = sum(!status$converged),
    improper_replicates = sum(status$improper, na.rm = TRUE),
    replicate_status = status
  )
}

# finally a wrapper for the bootstrap
bootstrap_sem <- function(df, spec) {

  # validate the arguments
  spec <- validate_sem_spec(spec)

  # check that this is a two-stage analysis
  if (!isTRUE(spec$bootstrap$requires_bootstrap)) {
    stop("bootstrap_sem() should only be used for two-stage analyses.")
  }

  # set the number of samples and the seed
  R <- spec$bootstrap$R
  set.seed(spec$bootstrap$seed)

  # fit once to get the paths
  original_fit <- fit_sem_once(df, spec)

  # vectors to store results
  rep_results <- vector("list", R)
  boot_draws <- vector("list", R)

  # loop over bootstrap replicates
  for (b in seq_len(R)) {
    boot_df <- draw_bootstrap_sample(df)

    fit_b <- tryCatch(
      fit_sem_once(boot_df, spec),
      error = function(e) e
    )

    if (inherits(fit_b, "error")) {
      rep_results[[b]] <- list(
        success = FALSE,
        converged = FALSE,
        improper = NA
      )
      boot_draws[[b]] <- NULL
      next
    }

    rep_results[[b]] <- list(
      success = TRUE,
      converged = isTRUE(fit_b$diagnostics$converged[1]),
      improper = isTRUE(fit_b$diagnostics$improper[1])
    )

    boot_draws[[b]] <- extract_bootstrap_targets(fit_b)
  }

  final_paths <- summarise_bootstrap_paths(
    original_paths = original_fit$paths,
    bootstrap_draws = boot_draws
  )

  failure_summary <- summarise_bootstrap_failures(rep_results, R)

  list(
    spec = spec,
    original_fit = original_fit,
    bootstrap_draws = boot_draws,
    final_paths = final_paths,
    failure_summary = failure_summary,
    diagnostics = original_fit$diagnostics,
    syntax = original_fit$syntax,
    preprocess_info = original_fit$preprocess_info
  )
}

# ------------------------------ 11) top level wrapper ------------------------------------------
# basically to decide which function to call
fit_sem_pipeline <- function(df, spec) {

  spec <- validate_sem_spec(spec)

  if (isTRUE(spec$bootstrap$requires_bootstrap)) {
    return(bootstrap_sem(df, spec))
  }

  fit <- fit_sem_once(df, spec)

  fit$paths$se_model <- fit$paths$se
  fit$paths$se_boot <- NA_real_
  fit$paths$se <- fit$paths$se_model

  fit
}

# --------------------------------- 12) write these into fitting functions ----------------------------
# CLPM
fit_clpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # residualiser confounders
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # CLPM direct confounders
  model_k = 0,
  model_exclude = NULL,
  model_confounder_order = 0,

  # bootstrap
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb fitting stage
  xgb_tuning = NULL,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning stage
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan
  estimator = "ML",
  fixed_x = FALSE
) {

  residualiser <- match.arg(residualiser)

  spec <- make_sem_spec(
    model_type = "clpm",
    T = T,
    residualiser = residualiser,

    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,

    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,

    model_k = model_k,
    model_exclude = model_exclude,
    model_confounder_order = model_confounder_order,

    free_loadings = FALSE,

    estimator = estimator,
    fixed_x = fixed_x,

    bootstrap_R = bootstrap_R,
    bootstrap_seed = bootstrap_seed,

    xgb_tuning = xgb_tuning,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_seed = xgb_seed,

    xgb_tuning_grid = xgb_tuning_grid,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds
  )

  fit_sem_pipeline(df, spec)
}

# RICLPM
fit_riclpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # residualiser confounders
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # measurement
  free_loadings = FALSE,

  # bootstrap
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb fitting stage
  xgb_tuning = NULL,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning stage
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan
  estimator = "ML",
  fixed_x = FALSE
) {

  residualiser <- match.arg(residualiser)

  spec <- make_sem_spec(
    model_type = "riclpm",
    T = T,
    residualiser = residualiser,

    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,

    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,

    model_k = 0,
    model_exclude = NULL,
    model_confounder_order = 0,

    free_loadings = free_loadings,

    estimator = estimator,
    fixed_x = fixed_x,

    bootstrap_R = bootstrap_R,
    bootstrap_seed = bootstrap_seed,

    xgb_tuning = xgb_tuning,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_seed = xgb_seed,

    xgb_tuning_grid = xgb_tuning_grid,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds
  )

  fit_sem_pipeline(df, spec)
}

# DPM
fit_dpm <- function(
  df,
  T,
  residualiser = c("none", "linear", "xgb"),

  x_prefix = "x",
  y_prefix = "y",
  c_prefix = "c",

  # residualiser confounders
  resid_k = NULL,
  resid_exclude = NULL,
  resid_interaction_order = 1,

  # measurement
  free_loadings = FALSE,

  # bootstrap
  bootstrap_R = NULL,
  bootstrap_seed = 123,

  # xgb fitting stage
  xgb_tuning = NULL,
  xgb_oof_folds = 2,
  xgb_nthread = 1,
  xgb_seed = 123,

  # xgb tuning stage
  xgb_tuning_grid = NULL,
  xgb_cv_folds = 5,
  xgb_nrounds_max = 400,
  xgb_early_stopping_rounds = 20,

  # lavaan
  estimator = "ML",
  fixed_x = FALSE
) {

  residualiser <- match.arg(residualiser)

  spec <- make_sem_spec(
    model_type = "dpm",
    T = T,
    residualiser = residualiser,

    x_prefix = x_prefix,
    y_prefix = y_prefix,
    c_prefix = c_prefix,

    resid_k = resid_k,
    resid_exclude = resid_exclude,
    resid_interaction_order = resid_interaction_order,

    model_k = 0,
    model_exclude = NULL,
    model_confounder_order = 0,

    free_loadings = free_loadings,

    estimator = estimator,
    fixed_x = fixed_x,

    bootstrap_R = bootstrap_R,
    bootstrap_seed = bootstrap_seed,

    xgb_tuning = xgb_tuning,
    xgb_oof_folds = xgb_oof_folds,
    xgb_nthread = xgb_nthread,
    xgb_seed = xgb_seed,

    xgb_tuning_grid = xgb_tuning_grid,
    xgb_cv_folds = xgb_cv_folds,
    xgb_nrounds_max = xgb_nrounds_max,
    xgb_early_stopping_rounds = xgb_early_stopping_rounds
  )

  fit_sem_pipeline(df, spec)
}