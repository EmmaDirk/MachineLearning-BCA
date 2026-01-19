# in a true application of XBG residualisation, one would want to tune the hyperparameters
# however, this is computationally very costly. 
# therefore, we provide a function that can be run once per study design to find the optimal hyperparameters
# and recycle these for all simulations in that design, recreating how tuning would happen in practice
# ------------------------------------------------------------------------------------------------------------

tune_xgb_once <- function(df,
                          k,
                          x_prefix = "x",
                          y_prefix = "y",
                          c_prefix = "c",
                          nfold = 5,
                          nrounds_max = 6000,                 # (CHANGED) allow more rounds for smaller eta
                          early_stopping_rounds = 100,        # (CHANGED) reduce "too-early stop"
                          max_grid = 600,                     # (NEW) slightly larger search budget
                          seed = 1) {                         # (NEW) reproducible grid + folds

  # convert to data frame
  df <- as.data.frame(df)

  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # only linear confounders are observed
  c_cols <- paste0(c_prefix, 1:k)

  # confounder matrix (standardised for XGB stability)
  C <- scale(as.matrix(df[c_cols]))

  # (NEW) helper: tune on a stacked target with grouped folds by person
  tune_one <- function(y_stack, C_stack, n_person) {

    dtrain <- xgboost::xgb.DMatrix(data = C_stack, label = y_stack)

    # (NEW) grouped folds (avoid leakage when the same person's C row is repeated across waves)
    person_id <- rep(seq_len(n_person), times = nrow(C_stack) / n_person)

    set.seed(seed)
    fold_id_person <- sample(rep(1:nfold, length.out = n_person))

    folds <- lapply(1:nfold, function(f) which(fold_id_person[person_id] == f))

    # grid of hyperparameters
    grid <- expand.grid(
      max_depth        = c(2:10),
      min_child_weight = c(1, 2, 5, 10, 20),
      eta              = c(0.003, 0.005, 0.01, 0.02, 0.05, 0.10),
      subsample        = c(0.5, 0.6, 0.8, 1.0),
      colsample_bytree = c(0.5, 0.6, 0.8, 1.0),
      gamma            = c(0, 0.1, 0.5, 1, 2, 5),
      lambda           = c(0.5, 1, 2, 5),
      alpha            = c(0, 0.001, 0.01, 0.1),
      KEEP.OUT.ATTRS   = FALSE
    )

    # trim the grid to keep tuning feasible
    if (nrow(grid) > max_grid) {
      set.seed(seed)
      grid <- grid[sample.int(nrow(grid), max_grid), , drop = FALSE]
    }

    best_rmse    <- Inf
    best_params  <- NULL
    best_nrounds <- NULL

    # create a progress bar for the CV grid
    # (this one shows elapsed time + estimated time left)
    pb <- pbapply::timerProgressBar(min = 0, max = nrow(grid), style = 3)

    for (i in seq_len(nrow(grid))) {

      # update progress bar
      pbapply::setTimerProgressBar(pb, i)

      params <- as.list(grid[i, ])

      cv <- xgboost::xgb.cv(
        data = dtrain,
        nrounds = nrounds_max,
        folds = folds,                          # (NEW) grouped folds to avoid leakage
        early_stopping_rounds = early_stopping_rounds,
        verbose = 0,
        params = c(
          list(
            objective = "reg:squarederror",
            eval_metric = "rmse",
            booster = "gbtree",
            tree_method = "hist",

            # let xgboost use multiple threads during tuning
            nthread = max(1, parallel::detectCores() - 1)
          ),
          params
        )
      )

      rmse <- cv$evaluation_log$test_rmse_mean[cv$best_iteration]

      if (rmse < best_rmse) {
        best_rmse    <- rmse
        best_params  <- params
        best_nrounds <- cv$best_iteration
      }
    }

    # close the progress bar
    close(pb)

    list(
      params  = best_params,
      nrounds = best_nrounds,
      rmse    = best_rmse
    )
  }

  # message so it is clear why the simulation progress bar has not started yet
  cat("\nTuning the XGB model (one-time CV)...\n")

  # (NEW) tune X and Y separately
  # stack all X waves to tune on more rows
  y_stack_x <- unlist(df[x_cols])
  C_stack_x <- C[rep(seq_len(nrow(C)), times = length(x_cols)), , drop = FALSE]

  cat("\n - Tuning X residualiser...\n")
  tuned_x <- tune_one(y_stack_x, C_stack_x, n_person = nrow(C))

  # stack all Y waves to tune on more rows
  y_stack_y <- unlist(df[y_cols])
  C_stack_y <- C[rep(seq_len(nrow(C)), times = length(y_cols)), , drop = FALSE]

  cat("\n - Tuning Y residualiser...\n")
  tuned_y <- tune_one(y_stack_y, C_stack_y, n_person = nrow(C))

  # newline so the next progress bar starts on a clean line
  cat("\n")

  list(
    params_X  = tuned_x$params,
    nrounds_X = tuned_x$nrounds,
    rmse_X    = tuned_x$rmse,

    params_Y  = tuned_y$params,
    nrounds_Y = tuned_y$nrounds,
    rmse_Y    = tuned_y$rmse
  )
}