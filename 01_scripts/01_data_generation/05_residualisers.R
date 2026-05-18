# =================================================================================================
#
# In this script we outline essentially what makes BCA SEM exactly that.
# We specify the functions that residualise variables X and Y with respect to confounders C.
#
# X ~ C + E, where E_x is the residual.
# Y ~ C + E, where E_y is the residual.
#
# We are interested in replacing X and Y with their residuals E.
# We can use multiple functions to model X ~ C, where some are more flexible than others.
# In this script we outline the following techniques to model that function:
# - 1) linear regression
# - 2) Extreme Gradient Boosting (XGBoost)
# - 3) Elastic Net via glmnet
# =================================================================================================

# ---- fold helpers --------------------------------------------------------------------------------

# Helper function to create grouped OOF folds.
# If duplicated bootstrap rows share the same original id, they are forced into the same fold.

make_group_folds <- function(id, oof_folds = 2, seed = 123) {

  # Stop if fold count is not valid.

  if (!is.numeric(oof_folds) || length(oof_folds) != 1 || is.na(oof_folds) || oof_folds < 2)
    stop("'oof_folds' must be a single integer >= 2.")

  # Coerce to integer.

  oof_folds <- as.integer(oof_folds)

  # Stop if id is missing.

  if (is.null(id))
    stop("'id' must not be NULL.")

  # Coerce to character to avoid factor / numeric indexing issues.

  id_chr <- as.character(id)

  # Get unique ids.

  u <- unique(id_chr)

  # Stop if there are fewer unique ids than folds.

  if (length(u) < oof_folds)
    stop("Need at least as many unique ids as OOF folds.")

  # Assign folds at the original-id level.

  set.seed(seed)
  fold_u <- sample(rep(seq_len(oof_folds), length.out = length(u)))
  names(fold_u) <- u

  # Map the original-id fold assignment back to all rows.

  as.integer(fold_u[id_chr])
}


# Helper to create ordinary OOF folds or grouped OOF folds for bootstrap resamples.

make_oof_fold_vector <- function(df, oof_folds = 2, seed = 123) {

  # Save sample size.

  n <- nrow(df)

  # If bootstrap preserved original row ids, keep duplicated copies in the same fold.

  if (".id_orig" %in% names(df)) {
    return(make_group_folds(df$.id_orig, oof_folds = oof_folds, seed = seed))
  }

  set.seed(seed)

  # Create folds.

  sample(rep(seq_len(oof_folds), length.out = n))
}


# ---- prediction metrics --------------------------------------------------------------------------

# Helper to compute prediction quality from out-of-fold predictions.
# We store both the MSE and the R^2 because they answer slightly different questions:
# - MSE is the average squared prediction error
# - R^2 is the fraction of variance in the target explained by the OOF predictions

compute_oof_metrics <- function(y, pred) {

  # Basic prediction error.

  mse <- mean((y - pred)^2, na.rm = TRUE)

  # Use the sample variance of the observed target.

  vy <- stats::var(y, na.rm = TRUE)

  # If the target has zero or undefined variance, R^2 is not defined.

  if (is.na(vy) || !is.finite(vy) || vy <= 0) {
    r2 <- NA_real_
  } else {
    r2 <- 1 - mse / vy
  }

  list(
    mse = as.numeric(mse),
    r2 = as.numeric(r2)
  )
}


# Helper to extract the wave number from names like x1, x2, y3, ...

extract_wave_number <- function(varname) {

  as.integer(sub("^[^0-9]+", "", varname))
}


# Helper to build one empty T-row ML metric frame.

make_panel_ml_metric_frame <- function(x_cols, y_cols) {

  # Prepare one T-row metric frame that stores the OOF prediction diagnostics
  # separately for X and Y at each observed wave.

  wave_x <- extract_wave_number(x_cols)
  wave_y <- extract_wave_number(y_cols)
  T_panel <- max(c(wave_x, wave_y), na.rm = TRUE)

  data.frame(
    T = seq_len(T_panel),
    mse_x = rep(NA_real_, T_panel),
    r2_x = rep(NA_real_, T_panel),
    mse_y = rep(NA_real_, T_panel),
    r2_y = rep(NA_real_, T_panel),
    stringsAsFactors = FALSE
  )
}


# ---- confounder-column selection -----------------------------------------------------------------

# Select observed confounder columns.
# This helper separates what the analyst gets to observe from how the residualiser models it.
#
# interaction_order determines which Appendix A feature blocks are visible:
# - 1 = base confounders only, such as C1 and C2
# - 2 = base confounders plus standardized two-way terms, such as Z1.2
# - 3 = base confounders plus standardized two-way and three-way terms, such as T1.2.3
#
# The data have already been made lavaan-safe before this helper is called, so interaction names
# usually contain dots rather than colons. The parser accepts both forms for robustness.

select_observed_confounder_columns <- function(df,
                                               k = NULL,
                                               c_prefix = "C",
                                               exclude = NULL,
                                               interaction_order = 1) {

  # Only main effects, two-way terms, and three-way terms are supported.

  if (!(interaction_order %in% c(1, 2, 3)))
    stop("'interaction_order' must be 1, 2, or 3.")

  # This parser extracts the integer indices from C, Z, or T feature names.

  parse_feature_name <- function(x) {
    parts <- strsplit(gsub(":", ".", x, fixed = TRUE), "\\.")[[1]]
    parts[1] <- sub("^[CcZzTt]", "", parts[1])
    as.integer(parts)
  }

  # This classifier maps a column name to an Appendix A feature block.

  feature_order <- function(x) {
    if (grepl(paste0("^", c_prefix, "\\d+$"), x)) return(1L)
    if (grepl("^Z\\d+([.:]\\d+)$", x)) return(2L)
    if (grepl("^T\\d+([.:]\\d+){2}$", x)) return(3L)
    NA_integer_
  }

  # A name is allowed when it belongs to a requested block and uses valid confounder indices.

  is_allowed_observed_name <- function(x) {

    order <- feature_order(x)

    if (is.na(order) || order > interaction_order) {
      return(FALSE)
    }

    idx <- parse_feature_name(x)

    if (any(is.na(idx))) {
      return(FALSE)
    }

    if (!is.null(k) && any(idx > k)) {
      return(FALSE)
    }

    TRUE
  }

  # Keep the selected confounder columns in the original data-column order.

  c_cols <- names(df)[vapply(names(df), is_allowed_observed_name, logical(1))]

  # Stop early if no observed confounder columns match the requested visibility level.

  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # Exclusions are applied after the visible confounder set has been selected.

  if (!is.null(exclude)) {

    if (!is.character(exclude))
      stop("'exclude' must be a character vector, e.g. exclude = c('C1', 'C2')")

    missing_exclude <- setdiff(exclude, names(df))

    if (length(missing_exclude) > 0)
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))

    c_cols <- setdiff(c_cols, exclude)
  }

  # At least one confounder must remain after exclusions.

  if (length(c_cols) == 0)
    stop("No confounders left after exclusion.")

  c_cols
}


# ---- design-matrix helpers -----------------------------------------------------------------------

# Helper to build the right-hand side of a linear-style formula.
# Once the observed confounder columns have been selected, the linear residualiser
# uses exactly those observed columns linearly. It does not generate extra terms.

build_linear_rhs <- function(c_cols) {

  paste(c_cols, collapse = " + ")
}


# Helper to build the confounder design matrix for xgboost.
# XGBoost sees exactly the observed confounder columns selected upstream.
# It does not receive an additional hand-built interaction expansion here.

build_xgb_confounder_matrix <- function(df, c_cols) {

  # Build design matrix with no intercept from the observed confounder columns only.

  X <- as.matrix(df[, c_cols, drop = FALSE])

  # Return numeric matrix.

  X
}


# Helper to build the elastic-net design matrix.
# Here we include:
# - linear terms
# - squared terms
# - cubic terms
# - all interaction terms up to third order across the observed confounder columns
#
# interaction_order does not control this modeling layer. It only controls which
# confounder columns are observed upstream. Once those columns are visible, Elastic Net
# always casts the same wide net over them.

build_enet_confounder_matrix <- function(df, c_cols) {

  # Build one cubic polynomial block for each observed confounder column.

  poly_terms <- vapply(
    c_cols,
    FUN = function(cc) sprintf("poly(%s, 3, raw = TRUE)", cc),
    FUN.VALUE = character(1)
  )

  # Build the full formula up to three-way interactions across those cubic blocks.

  rhs <- paste0("(", paste(poly_terms, collapse = " + "), ")^3")

  # Build design matrix with no intercept.

  X <- stats::model.matrix(stats::as.formula(paste("~", rhs, "- 1")), data = df)

  # Return numeric matrix.

  X
}


# ---- linear residualiser -------------------------------------------------------------------------

residualise_panel_linearC <- function(df,
                                      k = NULL,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "C",
                                      exclude = NULL,
                                      interaction_order = 1,
                                      oof_folds = 2,
                                      seed = 123) {

  # Convert to data frame.

  df <- as.data.frame(df)

  # Get column names for x and y variables.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Choose confounders.

  c_cols <- select_observed_confounder_columns(
    df = df,
    k = k,
    c_prefix = c_prefix,
    exclude = exclude,
    interaction_order = interaction_order
  )

  # Build the right-hand side of the formula.

  rhs <- build_linear_rhs(c_cols = c_cols)

  # Prepare one T-row metric frame that stores the OOF prediction diagnostics
  # separately for X and Y at each observed wave.

  ml_metrics <- make_panel_ml_metric_frame(x_cols = x_cols, y_cols = y_cols)

  # Create folds.

  folds <- make_oof_fold_vector(df = df, oof_folds = oof_folds, seed = seed)

  # Save sample size.

  n <- nrow(df)

  # Function to residualise a single variable against the chosen confounder formula.

  residualise_one <- function(varname) {

    # Create formula varname ~ confounders.

    fml <- stats::as.formula(paste(varname, "~", rhs))

    # Create empty vector for out-of-fold predictions.

    pred <- numeric(n)

    # Fit out-of-fold linear models.

    for (f in seq_len(oof_folds)) {

      # Grab training and test sets.

      train <- folds != f
      test  <- folds == f

      # Fit linear model on the training fold only.

      fit <- stats::lm(fml, data = df[train, , drop = FALSE])

      # Predict on the held-out fold.

      pred[test] <- stats::predict(fit, newdata = df[test, , drop = FALSE])
    }

    # Compute prediction diagnostics on the OOF predictions.

    metrics <- compute_oof_metrics(df[[varname]], pred)

    list(
      residual = df[[varname]] - pred,
      mse = metrics$mse,
      r2 = metrics$r2
    )
  }

  # Residualise every observed wave, including wave 1.

  for (x in x_cols) {

    out_x <- residualise_one(x)

    # Replace the x column with its residuals.

    df[[x]] <- out_x$residual

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(x)
    ml_metrics$mse_x[ml_metrics$T == wave] <- out_x$mse
    ml_metrics$r2_x[ml_metrics$T == wave] <- out_x$r2
  }

  # Same for y.

  for (y in y_cols) {

    out_y <- residualise_one(y)

    # Replace the y column with its residuals.

    df[[y]] <- out_y$residual

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(y)
    ml_metrics$mse_y[ml_metrics$T == wave] <- out_y$mse
    ml_metrics$r2_y[ml_metrics$T == wave] <- out_y$r2
  }

  # Return the residualised data frame together with the OOF metrics.

  list(
    data = df,
    ml_metrics = ml_metrics
  )
}


# ---- XGB tuning ----------------------------------------------------------------------------------

# We tune XGBoost twice for X and twice for Y:
# - early tuning: x2 and y2, used for waves 1, 2, and 3
# - late tuning:  x4 and y4, used for waves 4 and 5
#
# This is useful for stepwise scenarios where the confounding mechanism changes after wave 3.

tune_residualise_panel_xgb <- function(
  df,                                       # data frame
  k = NULL,                                 # number of confounders
  x_prefix = "x",                           # prefix for x variables
  y_prefix = "y",                           # prefix for y variables
  c_prefix = "C",                           # prefix for C variables
  exclude = NULL,                           # confounders to exclude
  interaction_order = 1,                    # controls which observed confounder columns are visible
  tuning_grid = NULL,                       # COST: grid of hyperparameters to try
  cv_folds = 5,                             # COST: number of CV folds
  nrounds_max = 400,                        # COST: maximum number of boosting iterations
  early_stopping_rounds = 20,               # COST: early stopping rounds
  early_tune_wave = 2L,                     # wave used for early X/Y tuning
  late_tune_wave = 4L,                      # wave used for late X/Y tuning
  late_start_wave = 4L,                     # first wave that uses late tuning
  nthread = 1,                              # number of threads for tuning
  seed = 123
){

  # Check that xgboost is installed.

  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("xgboost required.")

  # Convert to data frame.

  df <- as.data.frame(df)

  # Grab column names.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Normalize and check wave-control arguments.

  early_tune_wave <- as.integer(early_tune_wave[1])
  late_tune_wave  <- as.integer(late_tune_wave[1])
  late_start_wave <- as.integer(late_start_wave[1])

  if (is.na(early_tune_wave) || early_tune_wave < 1L) {
    stop("'early_tune_wave' must be a positive integer.")
  }

  if (is.na(late_tune_wave) || late_tune_wave < 1L) {
    stop("'late_tune_wave' must be a positive integer.")
  }

  if (is.na(late_start_wave) || late_start_wave < 1L) {
    stop("'late_start_wave' must be a positive integer.")
  }

  # Stop if the requested tuning waves are not present.

  x_early_var <- paste0(x_prefix, early_tune_wave)
  y_early_var <- paste0(y_prefix, early_tune_wave)
  x_late_var  <- paste0(x_prefix, late_tune_wave)
  y_late_var  <- paste0(y_prefix, late_tune_wave)

  required_tuning_vars <- c(x_early_var, y_early_var, x_late_var, y_late_var)
  missing_tuning_vars <- setdiff(required_tuning_vars, names(df))

  if (length(missing_tuning_vars) > 0) {
    stop(
      "The requested XGB tuning waves are not available. Missing variables: ",
      paste(missing_tuning_vars, collapse = ", ")
    )
  }

  # Choose confounders.

  c_cols <- select_observed_confounder_columns(
    df = df,
    k = k,
    c_prefix = c_prefix,
    exclude = exclude,
    interaction_order = interaction_order
  )

  # Build confounder matrix.

  X <- build_xgb_confounder_matrix(
    df = df,
    c_cols = c_cols
  )

  # Default grid.

  if (is.null(tuning_grid)) {
    tuning_grid <- expand.grid(
      eta = c(.05, .1),
      max_depth = c(2, 3, 4),
      min_child_weight = c(1, 5),
      subsample = c(.8, 1),
      colsample_bytree = c(.8, 1)
    )
  }

  # Helper to grab the best iteration robustly across xgboost versions.

  get_best_iteration_xgb_cv <- function(cv) {

    # Old-style location.

    iter <- cv$best_iteration

    # Fallback used by some newer versions.

    if (is.null(iter) || length(iter) == 0 || is.na(iter)) {
      if (!is.null(cv$early_stop) && !is.null(cv$early_stop$best_iteration)) {
        iter <- cv$early_stop$best_iteration
      }
    }

    # Final fallback: use the last available iteration.

    if (is.null(iter) || length(iter) == 0 || is.na(iter)) {
      iter <- nrow(cv$evaluation_log)
    }

    as.integer(iter[1])
  }

  # Helper to grab the best score robustly across xgboost versions.

  get_best_score_xgb_cv <- function(cv, iter) {

    # Available columns in the evaluation log.

    cols <- names(cv$evaluation_log)

    # Prefer the canonical column name first.

    score_col <- "test_rmse_mean"

    # Fallback to a more flexible search if needed.

    if (!(score_col %in% cols)) {
      score_col <- grep("^test.*rmse.*mean$", cols, value = TRUE)
    }

    if (length(score_col) == 0) {
      score_col <- grep("test.*rmse", cols, value = TRUE)
    }

    if (length(score_col) == 0) {
      stop("Could not find a test RMSE column in cv$evaluation_log.")
    }

    score <- cv$evaluation_log[[score_col[1]]][iter]

    if (length(score) == 0 || is.na(score) || !is.finite(score)) {
      stop("Invalid tuning score produced.")
    }

    as.numeric(score[1])
  }

  # Helper function.

  tune_target <- function(y, target_name = "target") {

    # Predict target from confounder matrix.

    dtrain <- xgboost::xgb.DMatrix(X, label = y)

    message(sprintf("Tuning XGBoost for %s", target_name))

    results <- pbapply::pblapply(
      X = seq_len(nrow(tuning_grid)),
      FUN = function(i) {

        params <- list(
          booster = "gbtree",
          objective = "reg:squarederror",
          eval_metric = "rmse",
          eta = tuning_grid$eta[i],
          max_depth = tuning_grid$max_depth[i],
          min_child_weight = tuning_grid$min_child_weight[i],
          subsample = tuning_grid$subsample[i],
          colsample_bytree = tuning_grid$colsample_bytree[i],
          nthread = nthread
        )

        set.seed(seed)

        # Run CV for each hyperparameter combination.

        out <- tryCatch({

          cv <- xgboost::xgb.cv(
            params = params,
            data = dtrain,
            nrounds = nrounds_max,
            nfold = cv_folds,
            early_stopping_rounds = early_stopping_rounds,
            verbose = 0
          )

          # Grab best iteration and score robustly.

          iter <- get_best_iteration_xgb_cv(cv)
          score <- get_best_score_xgb_cv(cv, iter)

          c(
            params,
            list(best_iter = iter, score = score)
          )

        }, error = function(e) {

          message(sprintf(
            "Tuning failed for %s at grid row %d: %s",
            target_name, i, conditionMessage(e)
          ))

          NULL
        })

        out
      }
    )

    # Keep only successful tuning results.

    results <- Filter(Negate(is.null), results)

    # Stop if every tuning attempt failed.

    if (length(results) == 0) {
      stop(sprintf("All XGBoost tuning runs failed for %s.", target_name))
    }

    # Pick best hyperparameter combination.

    scores <- sapply(results, function(x) x$score)
    valid_idx <- which(!is.na(scores) & is.finite(scores))

    if (length(valid_idx) == 0) {
      stop(sprintf("No valid XGBoost tuning scores were produced for %s.", target_name))
    }

    best_idx <- valid_idx[which.min(scores[valid_idx])]
    best <- results[[best_idx]]

    list(
      params = best[names(best) %in% c(
        "booster", "objective", "eval_metric", "eta", "max_depth",
        "min_child_weight", "subsample", "colsample_bytree", "nthread"
      )],
      nrounds = best$best_iter,
      score = best$score
    )
  }

  # Tuning targets.

  x_early_var <- paste0(x_prefix, early_tune_wave)
  y_early_var <- paste0(y_prefix, early_tune_wave)
  x_late_var  <- paste0(x_prefix, late_tune_wave)
  y_late_var  <- paste0(y_prefix, late_tune_wave)

  cat("\nTuning XGBoost using early and late wave targets...\n")
  cat(" - Early X tuning target:", x_early_var, "\n")
  cat(" - Early Y tuning target:", y_early_var, "\n")
  cat(" - Late  X tuning target:", x_late_var, "\n")
  cat(" - Late  Y tuning target:", y_late_var, "\n")
  cat(" - Late tuning starts at wave:", late_start_wave, "\n\n")

  tune_x_early <- tune_target(df[[x_early_var]], target_name = x_early_var)
  tune_y_early <- tune_target(df[[y_early_var]], target_name = y_early_var)
  tune_x_late  <- tune_target(df[[x_late_var]],  target_name = x_late_var)
  tune_y_late  <- tune_target(df[[y_late_var]],  target_name = y_late_var)

  # Store both the new wave-specific names and the old tune_x/tune_y names.
  # The old names preserve backward compatibility and point to the early tuning.

  list(
    confounders = c_cols,
    interaction_order = interaction_order,
    design_colnames = colnames(X),
    tune_x_early = tune_x_early,
    tune_y_early = tune_y_early,
    tune_x_late = tune_x_late,
    tune_y_late = tune_y_late,
    tune_x = tune_x_early,
    tune_y = tune_y_early,
    early_tune_wave = early_tune_wave,
    late_tune_wave = late_tune_wave,
    late_start_wave = late_start_wave,
    tuning_rule = sprintf(
      "waves 1-%d use wave-%d tuning; waves %d+ use wave-%d tuning",
      late_start_wave - 1L,
      early_tune_wave,
      late_start_wave,
      late_tune_wave
    )
  )
}


# ---- XGB residualiser ----------------------------------------------------------------------------

residualise_panel_xgb <- function(
  df,                                        # data frame
  tuning,                                    # tuning results
  k = NULL,                                  # number of confounders
  x_prefix = "x",                            # prefix for X variables
  y_prefix = "y",                            # prefix for Y variables
  c_prefix = "C",                            # prefix for C variables
  exclude = NULL,                            # confounders to exclude
  interaction_order = 1,                     # controls which observed confounder columns are visible
  oof_folds = 2,                             # number of OOF folds
  use_oof = TRUE,                            # whether to use out-of-fold predictions
  late_start_wave = NULL,                    # first wave that uses late tuning; NULL reads from tuning
  nthread = 1,                               # number of threads for fitting
  seed = 123                                 # random seed
){

  # Check that xgboost is installed.

  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("xgboost required.")

  # Convert to data frame.

  df <- as.data.frame(df)

  # Grab x and y variables.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Choose confounders.

  c_cols <- select_observed_confounder_columns(
    df = df,
    k = k,
    c_prefix = c_prefix,
    exclude = exclude,
    interaction_order = interaction_order
  )

  # Build confounder matrix.

  X <- build_xgb_confounder_matrix(
    df = df,
    c_cols = c_cols
  )

  # Save sample size.

  n <- nrow(X)

  # Prepare one T-row metric frame that stores the OOF prediction diagnostics
  # separately for X and Y at each observed wave.

  ml_metrics <- make_panel_ml_metric_frame(x_cols = x_cols, y_cols = y_cols)

  # Stop if OOF switch is not valid.

  if (!is.logical(use_oof) || length(use_oof) != 1 || is.na(use_oof))
    stop("'use_oof' must be TRUE or FALSE.")

  # Determine the wave split between early and late tuning.

  if (is.null(late_start_wave)) {
    late_start_wave <- tuning$late_start_wave
  }

  if (is.null(late_start_wave) || length(late_start_wave) == 0 || is.na(late_start_wave[1])) {
    late_start_wave <- 4L
  }

  late_start_wave <- as.integer(late_start_wave[1])

  if (is.na(late_start_wave) || late_start_wave < 1L) {
    stop("'late_start_wave' must be NULL or a positive integer.")
  }

  # Create folds.

  if (isTRUE(use_oof)) {
    folds <- make_oof_fold_vector(df = df, oof_folds = oof_folds, seed = seed)
  } else {
    folds <- NULL
  }

  # Train the model on the folds and return out-of-fold predictions.

  oof_predict <- function(y, params, nrounds) {

    # Enforce thread setting.

    params$nthread <- nthread

    # If requested, fit once on the full sample and predict the full sample.

    if (!isTRUE(use_oof)) {

      # Create DMatrix objects.

      dtrain <- xgboost::xgb.DMatrix(X, label = y)
      dtest  <- xgboost::xgb.DMatrix(X)

      # Train model.

      model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      # Grab predictions from the full sample.

      return(stats::predict(model, dtest))
    }

    # Create empty vector.

    pred <- numeric(n)

    # For each fold.

    for (f in seq_len(oof_folds)) {

      # Grab training and test sets.

      train <- folds != f
      test  <- folds == f

      # Create DMatrix objects.

      dtrain <- xgboost::xgb.DMatrix(X[train, , drop = FALSE], label = y[train])
      dtest  <- xgboost::xgb.DMatrix(X[test, , drop = FALSE])

      # Train model.

      model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      # Grab predictions from held-out fold.

      pred[test] <- stats::predict(model, dtest)
    }

    pred
  }

  # Determine which hyperparameters to use for each variable.
  # Waves before late_start_wave use early tuning; waves from late_start_wave onward use late tuning.

  get_spec <- function(var, prefix) {

    wave <- extract_wave_number(var)
    use_late <- !is.na(wave) && wave >= late_start_wave

    if (prefix == x_prefix) {

      if (isTRUE(use_late)) {
        if (!is.null(tuning$tune_x_late)) return(tuning$tune_x_late)
        if (!is.null(tuning$final$x_late)) return(tuning$final$x_late)
        stop(sprintf("No valid late X tuning specification found for wave >= %d.", late_start_wave))
      }

      if (!is.null(tuning$tune_x_early)) return(tuning$tune_x_early)
      if (!is.null(tuning$final$x_early)) return(tuning$final$x_early)
      if (!is.null(tuning$final$x_all)) return(tuning$final$x_all)
      if (!is.null(tuning$final$x_wave2plus)) return(tuning$final$x_wave2plus)
      if (!is.null(tuning$tune_x)) return(tuning$tune_x)
      stop("No valid early X tuning specification found.")

    } else {

      if (isTRUE(use_late)) {
        if (!is.null(tuning$tune_y_late)) return(tuning$tune_y_late)
        if (!is.null(tuning$final$y_late)) return(tuning$final$y_late)
        stop(sprintf("No valid late Y tuning specification found for wave >= %d.", late_start_wave))
      }

      if (!is.null(tuning$tune_y_early)) return(tuning$tune_y_early)
      if (!is.null(tuning$final$y_early)) return(tuning$final$y_early)
      if (!is.null(tuning$final$y_all)) return(tuning$final$y_all)
      if (!is.null(tuning$final$y_wave2plus)) return(tuning$final$y_wave2plus)
      if (!is.null(tuning$tune_y)) return(tuning$tune_y)
      stop("No valid early Y tuning specification found.")
    }
  }

  # Residualise X.

  for (x in x_cols) {

    # Get hyperparameter specification for this x variable.

    spec <- get_spec(x, x_prefix)

    # Get out-of-fold predictions from confounders.

    pred <- oof_predict(df[[x]], spec$params, spec$nrounds)

    # Compute prediction diagnostics on the OOF predictions.

    metrics <- compute_oof_metrics(df[[x]], pred)

    # Replace the x column with its residuals.

    df[[x]] <- df[[x]] - pred

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(x)
    ml_metrics$mse_x[ml_metrics$T == wave] <- metrics$mse
    ml_metrics$r2_x[ml_metrics$T == wave] <- metrics$r2
  }

  # Residualise Y.

  for (y in y_cols) {

    # Get hyperparameter specification for this y variable.

    spec <- get_spec(y, y_prefix)

    # Get out-of-fold predictions from confounders.

    pred <- oof_predict(df[[y]], spec$params, spec$nrounds)

    # Compute prediction diagnostics on the OOF predictions.

    metrics <- compute_oof_metrics(df[[y]], pred)

    # Replace the y column with its residuals.

    df[[y]] <- df[[y]] - pred

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(y)
    ml_metrics$mse_y[ml_metrics$T == wave] <- metrics$mse
    ml_metrics$r2_y[ml_metrics$T == wave] <- metrics$r2
  }

  # Return the residualised data frame together with the OOF metrics.

  list(
    data = df,
    ml_metrics = ml_metrics
  )
}


# ---- elastic-net tuning --------------------------------------------------------------------------

# We tune only once for X and once for Y, using wave 2 as the representative wave.
# The resulting X tuning is then reused for all X waves, and the resulting Y tuning
# is reused for all Y waves.

tune_residualise_panel_enet <- function(
  df,                                       # data frame
  k = NULL,                                 # number of confounders
  x_prefix = "x",                           # prefix for x variables
  y_prefix = "y",                           # prefix for y variables
  c_prefix = "C",                           # prefix for C variables
  exclude = NULL,                           # confounders to exclude
  interaction_order = 1,                    # interaction order
  alpha_grid = seq(0, 1, by = .1),          # COST: elastic-net mixing values to try
  cv_folds = 5,                             # COST: number of CV folds
  seed = 123
){

  # Check that glmnet is installed.

  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("glmnet required.")

  # Convert to data frame.

  df <- as.data.frame(df)

  # Grab column names.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Stop if not enough waves.

  if (length(x_cols) < 2 || length(y_cols) < 2)
    stop("Need at least 2 waves for tuning because tuning now uses wave 2.")

  # Choose confounders.

  c_cols <- select_observed_confounder_columns(
    df = df,
    k = k,
    c_prefix = c_prefix,
    exclude = exclude,
    interaction_order = interaction_order
  )

  # Build confounder matrix.

  X <- build_enet_confounder_matrix(
    df = df,
    c_cols = c_cols
  )

  # Helper function.

  tune_target <- function(y, target_name = "target") {

    message(sprintf("Tuning Elastic Net for %s", target_name))

    results <- pbapply::pblapply(
      X = seq_along(alpha_grid),
      FUN = function(i) {

        alpha_i <- alpha_grid[i]

        set.seed(seed)

        out <- tryCatch({

          cv_fit <- glmnet::cv.glmnet(
            x = X,
            y = y,
            family = "gaussian",
            alpha = alpha_i,
            nfolds = cv_folds,
            standardize = TRUE,
            intercept = TRUE,
            type.measure = "mse"
          )

          list(
            alpha = alpha_i,
            lambda = cv_fit$lambda.min,
            score = min(cv_fit$cvm, na.rm = TRUE)
          )

        }, error = function(e) {

          message(sprintf(
            "Tuning failed for %s at alpha %.3f: %s",
            target_name, alpha_i, conditionMessage(e)
          ))

          NULL
        })

        out
      }
    )

    # Keep only successful tuning results.

    results <- Filter(Negate(is.null), results)

    # Stop if every tuning attempt failed.

    if (length(results) == 0) {
      stop(sprintf("All Elastic Net tuning runs failed for %s.", target_name))
    }

    # Pick best hyperparameter combination.

    scores <- sapply(results, function(x) x$score)
    valid_idx <- which(!is.na(scores) & is.finite(scores))

    if (length(valid_idx) == 0) {
      stop(sprintf("No valid Elastic Net tuning scores were produced for %s.", target_name))
    }

    best_idx <- valid_idx[which.min(scores[valid_idx])]
    best <- results[[best_idx]]

    list(
      alpha = as.numeric(best$alpha),
      lambda = as.numeric(best$lambda),
      score = as.numeric(best$score)
    )
  }

  # Tune once for x2 and once for y2.

  x_tune_var <- paste0(x_prefix, 2)
  y_tune_var <- paste0(y_prefix, 2)

  cat("\nTuning Elastic Net using representative wave 2 targets...\n")
  cat(" - X tuning target:", x_tune_var, "\n")
  cat(" - Y tuning target:", y_tune_var, "\n\n")

  tune_x <- tune_target(df[[x_tune_var]], target_name = x_tune_var)
  tune_y <- tune_target(df[[y_tune_var]], target_name = y_tune_var)

  # Store the compact tuning result.
  # We now use one X tuning object and one Y tuning object for all waves.

  list(
    confounders = c_cols,
    interaction_order = interaction_order,
    design_colnames = colnames(X),
    tune_x = tune_x,
    tune_y = tune_y
  )
}


# ---- elastic-net residualiser --------------------------------------------------------------------

residualise_panel_enet <- function(
  df,                                        # data frame
  tuning,                                    # tuning results
  k = NULL,                                  # number of confounders
  x_prefix = "x",                            # prefix for X variables
  y_prefix = "y",                            # prefix for Y variables
  c_prefix = "C",                            # prefix for C variables
  exclude = NULL,                            # confounders to exclude
  interaction_order = 1,                     # controls which observed confounder columns are visible
  oof_folds = 2,                             # number of OOF folds
  seed = 123                                 # random seed
){

  # Check that glmnet is installed.

  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("glmnet required.")

  # Convert to data frame.

  df <- as.data.frame(df)

  # Grab x and y variables.

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # Choose confounders.

  c_cols <- select_observed_confounder_columns(
    df = df,
    k = k,
    c_prefix = c_prefix,
    exclude = exclude,
    interaction_order = interaction_order
  )

  # Build confounder matrix.

  X <- build_enet_confounder_matrix(
    df = df,
    c_cols = c_cols
  )

  # Save sample size.

  n <- nrow(X)

  # Prepare one T-row metric frame that stores the OOF prediction diagnostics
  # separately for X and Y at each observed wave.

  ml_metrics <- make_panel_ml_metric_frame(x_cols = x_cols, y_cols = y_cols)

  # Create folds.

  folds <- make_oof_fold_vector(df = df, oof_folds = oof_folds, seed = seed)

  # Train the model on the folds and return out-of-fold predictions.

  oof_predict <- function(y, alpha, lambda) {

    # Create empty vector.

    pred <- numeric(n)

    # For each fold.

    for (f in seq_len(oof_folds)) {

      # Grab training and test sets.

      train <- folds != f
      test  <- folds == f

      # Fit model on the training fold only.

      fit <- glmnet::glmnet(
        x = X[train, , drop = FALSE],
        y = y[train],
        family = "gaussian",
        alpha = alpha,
        lambda = lambda,
        standardize = TRUE,
        intercept = TRUE
      )

      # Grab predictions from held-out fold.

      pred[test] <- as.numeric(stats::predict(
        fit,
        newx = X[test, , drop = FALSE],
        s = lambda
      ))
    }

    pred
  }

  # Determine which hyperparameters to use for each variable.

  get_spec <- function(prefix) {

    if (prefix == x_prefix) {

      if (!is.null(tuning$tune_x)) return(tuning$tune_x)
      stop("No valid Elastic Net tuning specification found for X.")

    } else {

      if (!is.null(tuning$tune_y)) return(tuning$tune_y)
      stop("No valid Elastic Net tuning specification found for Y.")
    }
  }

  # Residualise X.

  for (x in x_cols) {

    # Get hyperparameter specification for this x variable.

    spec <- get_spec(x_prefix)

    # Get out-of-fold predictions from confounders.

    pred <- oof_predict(df[[x]], spec$alpha, spec$lambda)

    # Compute prediction diagnostics on the OOF predictions.

    metrics <- compute_oof_metrics(df[[x]], pred)

    # Replace the x column with its residuals.

    df[[x]] <- df[[x]] - pred

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(x)
    ml_metrics$mse_x[ml_metrics$T == wave] <- metrics$mse
    ml_metrics$r2_x[ml_metrics$T == wave] <- metrics$r2
  }

  # Residualise Y.

  for (y in y_cols) {

    # Get hyperparameter specification for this y variable.

    spec <- get_spec(y_prefix)

    # Get out-of-fold predictions from confounders.

    pred <- oof_predict(df[[y]], spec$alpha, spec$lambda)

    # Compute prediction diagnostics on the OOF predictions.

    metrics <- compute_oof_metrics(df[[y]], pred)

    # Replace the y column with its residuals.

    df[[y]] <- df[[y]] - pred

    # Store the OOF metrics at the matching wave.

    wave <- extract_wave_number(y)
    ml_metrics$mse_y[ml_metrics$T == wave] <- metrics$mse
    ml_metrics$r2_y[ml_metrics$T == wave] <- metrics$r2
  }

  # Return the residualised data frame together with the OOF metrics.

  list(
    data = df,
    ml_metrics = ml_metrics
  )
}