# In this script we outline essentially what makes BCA SEM exactly that.
# We specify the functions that residualise variables X and Y with respect to confounders C.
#
# X ~ C + E, where E_x is the residual
# Y ~ C + E, where E_y is the residual
#
# We are interested in replacing X and Y with their residuals E.
# We can use multiple functions to model X ~ C, where some are more flexible than others.
# In this script we outline the following techniques to model that function.
# - 1) linear regression.
# - 2) Extreme Gradient Boosting (XGBoost).
# - TBA) Penalized Regression (Lasso).
# - TBA) Tabular Prior Data-Fitted Network (TabPFN).
# ------------------------------------------------------------------------------------------

# helper function to create grouped OOF folds
# if duplicated bootstrap rows share the same original id, they are forced into the same fold
make_group_folds <- function(id, oof_folds = 2, seed = 123) {

  # stop if fold count is not valid
  if (!is.numeric(oof_folds) || length(oof_folds) != 1 || is.na(oof_folds) || oof_folds < 2)
    stop("'oof_folds' must be a single integer >= 2.")

  # coerce to integer
  oof_folds <- as.integer(oof_folds)

  # stop if id is missing
  if (is.null(id))
    stop("'id' must not be NULL.")

  # coerce to character to avoid factor / numeric indexing issues
  id_chr <- as.character(id)

  # get unique ids
  u <- unique(id_chr)

  # stop if there are fewer unique ids than folds
  if (length(u) < oof_folds)
    stop("Need at least as many unique ids as OOF folds.")

  # assign folds at the original-id level
  set.seed(seed)
  fold_u <- sample(rep(seq_len(oof_folds), length.out = length(u)))
  names(fold_u) <- u

  # map the original-id fold assignment back to all rows
  as.integer(fold_u[id_chr])
}

# ----------------------------------- 1) linear residualiser -------------------------------
residualise_panel_linearC <- function(df,
                                      k = NULL,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c",
                                      exclude = NULL,
                                      interaction_order = 1,
                                      oof_folds = 2,
                                      seed = 123) {

  # convert to data frame
  df <- as.data.frame(df)

  # get column names for x and y variables
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # choose confounders
  if (is.null(k)) {

    # use all confounder columns that match the prefix
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)

  } else {

    # use only confounders c1 to ck
    c_cols <- paste0(c_prefix, 1:k)
  }

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0)
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))

  # exclude selected confounders if requested
  if (!is.null(exclude)) {

    # stop if exclude is not given as character vector of column names
    if (!is.character(exclude))
      stop("'exclude' must be a character vector, e.g. exclude = c('c1', 'c2')")

    # stop if excluded confounders are not found in the data
    missing_exclude <- setdiff(exclude, names(df))
    if (length(missing_exclude) > 0)
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))

    # remove excluded confounders from the confounder set
    c_cols <- setdiff(c_cols, exclude)
  }

  # stop if no confounders remain after exclusion
  if (length(c_cols) == 0)
    stop("No confounders left after exclusion.")

  # stop if interaction order is not valid
  if (!(interaction_order %in% c(1, 2, 3)))
    stop("'interaction_order' must be 1, 2, or 3.")

  # build the right-hand side of the formula
  if (interaction_order == 1) {

    # main effects only: c1 + c2 + c3 + ...
    rhs <- paste(c_cols, collapse = " + ")

  } else {

    # include all effects up to the chosen interaction order
    # e.g. (c1 + c2 + c3)^2 gives main effects + all 2-way interactions
    # e.g. (c1 + c2 + c3)^3 gives main effects + all 2-way + all 3-way interactions
    rhs <- paste0("(", paste(c_cols, collapse = " + "), ")^", interaction_order)
  }

  # save sample size
  n <- nrow(df)

  # -------- create folds --------

  # if bootstrap preserved original row ids, keep duplicated copies in the same fold
  if (".id_orig" %in% names(df)) {

    # create grouped folds at the original-id level
    folds <- make_group_folds(df$.id_orig, oof_folds = oof_folds, seed = seed)

  } else {

    set.seed(seed)

    # create folds
    folds <- sample(rep(1:oof_folds, length.out = n))
  }

  # function to residualise a single variable against the chosen confounder formula
  residualise_one <- function(varname) {

    # create formula varname ~ confounders
    fml <- as.formula(paste(varname, "~", rhs))

    # create empty vector for out-of-fold predictions
    pred <- numeric(n)

    # fit out-of-fold linear models
    for (f in 1:oof_folds) {

      # grab training and test sets
      train <- folds != f
      test  <- folds == f

      # fit linear model on the training fold only
      fit <- lm(fml, data = df[train, , drop = FALSE])

      # predict on the held-out fold
      pred[test] <- predict(fit, newdata = df[test, , drop = FALSE])
    }

    # return OOF residuals
    df[[varname]] - pred
  }

  # IMPORTANT:
  # residualise every observed wave, including wave 1
  for (x in x_cols) {
    df[[x]] <- residualise_one(x)
  }

  # same for y
  for (y in y_cols) {
    df[[y]] <- residualise_one(y)
  }

  # return the residualised data frame
  df
}

# ----------------------------------- 2) xgb residualiser -------------------------------

# helper function to build the confounder design matrix for xgboost
# this is the key fix:
# - interaction_order = 1 gives only the main effects
# - interaction_order = 2 gives main effects + all 2-way interactions
# - interaction_order = 3 gives main effects + all 2-way + all 3-way interactions

build_xgb_confounder_matrix <- function(df, c_cols, interaction_order = 1) {

  # stop if interaction order is not valid
  if (!(interaction_order %in% c(1, 2, 3)))
    stop("'interaction_order' must be 1, 2, or 3.")

  # build the right-hand side of the formula
  if (interaction_order == 1) {

    # main effects only
    rhs <- paste(c_cols, collapse = " + ")

  } else {

    # include all effects up to the chosen interaction order
    rhs <- paste0("(", paste(c_cols, collapse = " + "), ")^", interaction_order)
  }

  # build design matrix with no intercept
  # example:
  # interaction_order = 1 --> c1, c2, c3
  # interaction_order = 2 --> c1, c2, c3, c1:c2, c1:c3, c2:c3
  # interaction_order = 3 --> main effects + 2-way + 3-way interactions
  X <- model.matrix(as.formula(paste("~", rhs, "- 1")), data = df)

  # return numeric matrix
  X
}

# ---------------------------- 2.1) tuning function ------------------------------------
# IMPORTANT:
# We now tune only once for X and once for Y, using wave 2 as the representative wave.
# The resulting X tuning is then reused for all X waves, and the resulting Y tuning
# is reused for all Y waves.

tune_residualise_panel_xgb <- function(
  df,                                       # data frame
  k = NULL,                                 # number of confounders
  x_prefix = "x",                           # prefix for x variables
  y_prefix = "y",                           # prefix for y variables
  c_prefix = "c",                           # prefix for c variables
  exclude = NULL,                           # confounders to exclude
  interaction_order = 1,                    # interaction order (so 1 uses only main effects)
  tuning_grid = NULL,                       # COST: grid of hyperparameters to try
  cv_folds = 5,                             # COST: number of CV folds
  nrounds_max = 400,                        # COST: maximum number of boosting iterations
  early_stopping_rounds = 20,               # COST: early stopping rounds
  nthread = 1,                              # number of threads for tuning
  seed = 123
){

  # check that xgboost is installed
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("xgboost required.")

  # convert to data frame
  df <- as.data.frame(df)

  # grab column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # stop if not enough waves
  if (length(x_cols) < 2 || length(y_cols) < 2)
    stop("Need at least 2 waves for tuning because tuning now uses wave 2.")

  # ---------------- confounder selection ----------------

  # choose confounders
  if (is.null(k)) {

    # use all confounder columns that match the prefix
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)

  } else {

    # use only confounders c1 to ck
    c_cols <- paste0(c_prefix, 1:k)
  }

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0)
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))

  # exclude selected confounders if requested
  if (!is.null(exclude)) {

    # stop if exclude is not given as character vector of column names
    if (!is.character(exclude))
      stop("'exclude' must be a character vector, e.g. exclude = c('c1', 'c2')")

    # stop if excluded confounders are not found in the data
    missing_exclude <- setdiff(exclude, names(df))
    if (length(missing_exclude) > 0)
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))

    # remove excluded confounders from the confounder set
    c_cols <- setdiff(c_cols, exclude)
  }

  # stop if no confounders remain after exclusion
  if (length(c_cols) == 0)
    stop("No confounders left after exclusion.")

  # build confounder matrix
  X <- build_xgb_confounder_matrix(
    df = df,
    c_cols = c_cols,
    interaction_order = interaction_order
  )

  # default grid
  if (is.null(tuning_grid)) {
    tuning_grid <- expand.grid(
      eta = c(.05, .1),
      max_depth = c(2, 3, 4),
      min_child_weight = c(1, 5),
      subsample = c(.8, 1),
      colsample_bytree = c(.8, 1)
    )
  }

  # helper to grab the best iteration robustly across xgboost versions
  get_best_iteration_xgb_cv <- function(cv) {

    # old-style location
    iter <- cv$best_iteration

    # fallback used by some newer versions
    if (is.null(iter) || length(iter) == 0 || is.na(iter)) {
      if (!is.null(cv$early_stop) && !is.null(cv$early_stop$best_iteration)) {
        iter <- cv$early_stop$best_iteration
      }
    }

    # final fallback: use the last available iteration
    if (is.null(iter) || length(iter) == 0 || is.na(iter)) {
      iter <- nrow(cv$evaluation_log)
    }

    as.integer(iter[1])
  }

  # helper to grab the best score robustly across xgboost versions
  get_best_score_xgb_cv <- function(cv, iter) {

    # available columns in the evaluation log
    cols <- names(cv$evaluation_log)

    # prefer the canonical column name first
    score_col <- "test_rmse_mean"

    # fallback to a more flexible search if needed
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

  # helper function
  tune_target <- function(y, target_name = "target") {

    # predict target from confounder matrix
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

        # run CV for each hyperparameter combination
        out <- tryCatch({

          cv <- xgboost::xgb.cv(
            params = params,
            data = dtrain,
            nrounds = nrounds_max,
            nfold = cv_folds,
            early_stopping_rounds = early_stopping_rounds,
            verbose = 0
          )

          # grab best iteration and score robustly
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

    # keep only successful tuning results
    results <- Filter(Negate(is.null), results)

    # stop if every tuning attempt failed
    if (length(results) == 0) {
      stop(sprintf("All XGBoost tuning runs failed for %s.", target_name))
    }

    # pick best hyperparameter combination
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

  # tune once for x2 and once for y2
  x_tune_var <- paste0(x_prefix, 2)
  y_tune_var <- paste0(y_prefix, 2)

  cat("\nTuning XGBoost using representative wave 2 targets...\n")
  cat(" - X tuning target:", x_tune_var, "\n")
  cat(" - Y tuning target:", y_tune_var, "\n\n")

  tune_x <- tune_target(df[[x_tune_var]], target_name = x_tune_var)
  tune_y <- tune_target(df[[y_tune_var]], target_name = y_tune_var)

  # store compact results, and also keep a backward-compatible 'final' block
  list(
    confounders = c_cols,
    interaction_order = interaction_order,
    design_colnames = colnames(X),

    tune_x = tune_x,
    tune_y = tune_y,

    final = list(
      x_all = tune_x,
      y_all = tune_y,

      # backward-compatible aliases
      x_wave1 = tune_x,
      x_wave2plus = tune_x,
      y_wave1 = tune_y,
      y_wave2plus = tune_y
    )
  )
}

# ---------------------- 2.2) training and residualising function ------------------------
residualise_panel_xgb <- function(
  df,                                        # data frame
  tuning,                                    # tuning results
  k = NULL,                                  # number of confounders
  x_prefix = "x",                            # prefix for X variables
  y_prefix = "y",                            # prefix for Y variables
  c_prefix = "c",                            # prefix for C variables
  exclude = NULL,                            # confounders to exclude
  interaction_order = 1,                     # interaction order
  oof_folds = 2,                             # number of OOF folds
  nthread = 1,                               # number of threads for fitting
  seed = 123                                 # random seed
){

  # check that xgboost is installed
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("xgboost required.")

  # convert to data frame
  df <- as.data.frame(df)

  # grab x and y variables
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # choose confounders
  if (is.null(k)) {

    # use all confounder columns that match the prefix
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)

  } else {

    # use only confounders c1 to ck
    c_cols <- paste0(c_prefix, 1:k)
  }

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0)
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))

  # exclude selected confounders if requested
  if (!is.null(exclude)) {

    # stop if exclude is not given as character vector of column names
    if (!is.character(exclude))
      stop("'exclude' must be a character vector, e.g. exclude = c('c1', 'c2')")

    # stop if excluded confounders are not found in the data
    missing_exclude <- setdiff(exclude, names(df))
    if (length(missing_exclude) > 0)
      stop("Excluded confounder columns not found: ", paste(missing_exclude, collapse = ", "))

    # remove excluded confounders from the confounder set
    c_cols <- setdiff(c_cols, exclude)
  }

  # stop if no confounders remain after exclusion
  if (length(c_cols) == 0)
    stop("No confounders left after exclusion.")

  # build confounder matrix
  X <- build_xgb_confounder_matrix(
    df = df,
    c_cols = c_cols,
    interaction_order = interaction_order
  )

  # save sample size
  n <- nrow(X)

  # -------- create folds --------

  # if bootstrap preserved original row ids, keep duplicated copies in the same fold
  if (".id_orig" %in% names(df)) {

    # create grouped folds at the original-id level
    folds <- make_group_folds(df$.id_orig, oof_folds = oof_folds, seed = seed)

  } else {

    set.seed(seed)

    # create folds
    folds <- sample(rep(1:oof_folds, length.out = n))
  }

  # train the model on the folds and return out-of-fold predictions
  oof_predict <- function(y, params, nrounds) {

    # create empty vector
    pred <- numeric(n)

    # for each fold
    for (f in 1:oof_folds) {

      # grab training and test sets
      train <- folds != f
      test  <- folds == f

      # create DMatrix objects
      dtrain <- xgboost::xgb.DMatrix(X[train, , drop = FALSE], label = y[train])
      dtest  <- xgboost::xgb.DMatrix(X[test, , drop = FALSE])

      # enforce thread setting
      params$nthread <- nthread

      # train model
      model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = nrounds,
        verbose = 0
      )

      # grab predictions from held-out fold
      pred[test] <- predict(model, dtest)
    }

    pred
  }

  # determine which hyperparameters to use for each variable
  # we now use one X tuning object for all X waves and one Y tuning object for all Y waves
  get_spec <- function(var, prefix) {

    if (prefix == x_prefix) {

      if (!is.null(tuning$final$x_all)) return(tuning$final$x_all)
      if (!is.null(tuning$final$x_wave2plus)) return(tuning$final$x_wave2plus)
      stop("No valid X tuning specification found.")

    } else {

      if (!is.null(tuning$final$y_all)) return(tuning$final$y_all)
      if (!is.null(tuning$final$y_wave2plus)) return(tuning$final$y_wave2plus)
      stop("No valid Y tuning specification found.")
    }
  }

  # residualise X
  for (x in x_cols) {

    # get hyperparameter specification for this x variable
    spec <- get_spec(x, x_prefix)

    # get out-of-fold predictions from confounders
    pred <- oof_predict(df[[x]], spec$params, spec$nrounds)

    # replace the x column with its residuals
    df[[x]] <- df[[x]] - pred
  }

  # residualise Y
  for (y in y_cols) {

    # get hyperparameter specification for this y variable
    spec <- get_spec(y, y_prefix)

    # get out-of-fold predictions from confounders
    pred <- oof_predict(df[[y]], spec$params, spec$nrounds)

    # replace the y column with its residuals
    df[[y]] <- df[[y]] - pred
  }

  # return the residualised data frame
  df
}