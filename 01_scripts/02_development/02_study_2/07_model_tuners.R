# In a true application of XGB residualisation, one would want to tune the hyperparameters
# However, it is computationally very constly to tune every simulation round. 
# Therefore, we provide a function that can be run once per study design to find the optimal hyperparameters
# and recycle these for all simulations in that design, recreating how tuning would happen in practice.
#
# This function supports tune profiles ("quick", "medium", "full") to control how expensive tuning is.
# - "quick": small grid + simpler trees (meant for testing your pipeline)
# - "medium": moderate grid + moderate complexity
# - "full": large grid + richer complexity (closest to "real" tuning)
# ------------------------------------------------------------------------------------------------------------

tune_xgb_once <- function(df,                                             # data frame with panel data
                          k,                                              # number of confounders
                          x_prefix = "x",                                 # prefix for X variables
                          y_prefix = "y",                                 # prefix for Y variables
                          c_prefix = "c",                                 # prefix for confounder variables
                          tune_profile = c("quick", "medium", "full"),    # tune profile (controls grid size and complexity)
                          nfold = NULL,                                   # number of CV folds (ovveriding profile default)
                          nrounds_max = NULL,                             # max number of boosting rounds (overriding profile default)
                          early_stopping_rounds = NULL,                   # early stopping rounds (overriding profile default)
                          max_grid = NULL,                                # max grid size (overriding profile default)
                          seed = 1) {                                     # random seed

  # match tune profile argument
  tune_profile <- match.arg(tune_profile)

  # profile defaults (that can be overridden)
  profile_defaults <- switch(

    # tune profiles
    tune_profile,

    # quick profile
    quick = list(

      # settings
      nfold = 3,                                        # 3-fold CV
      nrounds_max = 800,                                # max 800 rounds
      early_stopping_rounds = 40,                       # early stop after 40 rounds
      max_grid = 40,                                    # max 40 grid points
      grid = list(                                      # hyperparameter grid
        max_depth        = 2:4,                         # shallow trees
        min_child_weight = c(5, 10, 20),                # higher min child weight
        eta              = c(0.05, 0.10, 0.20),         # higher learning rate
        subsample        = c(0.7, 1.0),                 # fewer subsample options
        colsample_bytree = c(0.7, 1.0),                 # fewer colsample options
        gamma            = c(0, 1),                     # fewer gamma options
        lambda           = c(1, 2),                     # fewer lambda options
        alpha            = c(0, 0.01)                   # fewer alpha options
      )
    ),

    # medium profile
    medium = list(

      # settings
      nfold = 5,                                       # 5-fold CV
      nrounds_max = 2500,                              # max 2500 rounds
      early_stopping_rounds = 80,                      # early stop after 80 rounds
      max_grid = 200,                                  # max 200 grid points
      grid = list(                                     # hyperparameter grid
        max_depth        = 2:8,                        # deeper trees
        min_child_weight = c(1, 2, 5, 10, 20),         # moderate min child weight options
        eta              = c(0.01, 0.02, 0.05, 0.10),  # moderate learning rate options
        subsample        = c(0.6, 0.8, 1.0),           # moderate subsample options
        colsample_bytree = c(0.6, 0.8, 1.0),           # moderate colsample options
        gamma            = c(0, 0.1, 0.5, 1, 2),       # moderate gamma options
        lambda           = c(0.5, 1, 2, 5),            # moderate lambda options
        alpha            = c(0, 0.001, 0.01, 0.1)      # moderate alpha options
      )
    ),

    # full profile
    full = list(

      # settings
      nfold = 5,                                       # 5-fold CV
      nrounds_max = 6000,                              # max 6000 rounds
      early_stopping_rounds = 100,                     # early stop after 100 rounds
      max_grid = 600,                                  # max 600 grid points
      grid = list(                                     # hyperparameter grid
        max_depth        = 2:10,                       # full depth range
        min_child_weight = c(1, 2, 5, 10, 20),         # full min child weight options
        eta              = c(0.003, 0.005, 0.01, 0.02, # full learning rate options
          0.05, 0.10), 
        subsample        = c(0.5, 0.6, 0.8, 1.0),      # full subsample options
        colsample_bytree = c(0.5, 0.6, 0.8, 1.0),      # full colsample options
        gamma            = c(0, 0.1, 0.5, 1, 2, 5),    # full gamma options
        lambda           = c(0.5, 1, 2, 5),            # full lambda options
        alpha            = c(0, 0.001, 0.01, 0.1)      # full alpha options
      )
    )
  )

  # apply overrides if provided
  if (is.null(nfold)) nfold <- profile_defaults$nfold
  if (is.null(nrounds_max)) nrounds_max <- profile_defaults$nrounds_max
  if (is.null(early_stopping_rounds)) early_stopping_rounds <- profile_defaults$early_stopping_rounds
  if (is.null(max_grid)) max_grid <- profile_defaults$max_grid

  # convert to data frame
  df <- as.data.frame(df)

  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # only linear confounders are observed
  c_cols <- paste0(c_prefix, 1:k)

  # confounder matrix (standardised for XGB stability)
  C <- scale(as.matrix(df[c_cols]))

  # helper: tune on a stacked target with grouped folds by person
  tune_one <- function(y_stack, C_stack, n_person) {

    # create DMatrix
    dtrain <- xgboost::xgb.DMatrix(data = C_stack, label = y_stack)

    # grouped folds (avoid leakage when the same person's C row is repeated across waves)
    person_id <- rep(seq_len(n_person), times = nrow(C_stack) / n_person)

    # grouped CV folds
    set.seed(seed)

    # assign persons to folds
    fold_id_person <- sample(rep(1:nfold, length.out = n_person))

    # create folds based on person assignment
    folds <- lapply(1:nfold, function(f) which(fold_id_person[person_id] == f))

    # grid of hyperparameters (profile-dependent)
    grid <- expand.grid(
      max_depth        = profile_defaults$grid$max_depth,
      min_child_weight = profile_defaults$grid$min_child_weight,
      eta              = profile_defaults$grid$eta,
      subsample        = profile_defaults$grid$subsample,
      colsample_bytree = profile_defaults$grid$colsample_bytree,
      gamma            = profile_defaults$grid$gamma,
      lambda           = profile_defaults$grid$lambda,
      alpha            = profile_defaults$grid$alpha,
      KEEP.OUT.ATTRS   = FALSE
    )

    # trim the grid to keep tuning feasible
    # if the grid is larger than max_grid, sample without replacement
    if (nrow(grid) > max_grid) {

      # set random seed for reproducibility
      set.seed(seed)

      # sample without replacement
      grid <- grid[sample.int(nrow(grid), max_grid), , drop = FALSE]
    }

    # initialize best tracking variables
    best_rmse    <- Inf
    best_params  <- NULL
    best_nrounds <- NULL

    # progress bar for the CV grid
    pb <- pbapply::timerProgressBar(min = 0, max = nrow(grid), style = 3)
    
    # loop over grid rows
    for (i in seq_len(nrow(grid))) {
      
      # update progress bar
      pbapply::setTimerProgressBar(pb, i)
      
      # set hyperparameters
      params <- as.list(grid[i, ])

      # run CV
      cv <- xgboost::xgb.cv(

        # data and parameters
        data = dtrain,
        nrounds = nrounds_max,
        folds = folds,
        early_stopping_rounds = early_stopping_rounds,

        # silent
        verbose = 0,

        # fixed parameters
        params = c(
          list(
            objective = "reg:squarederror",
            eval_metric = "rmse",
            booster = "gbtree",
            tree_method = "hist",
            nthread = max(1, parallel::detectCores() - 1)
          ),

          # return the current grid row
          params
        )
      )

      # get best RMSE
      rmse <- cv$evaluation_log$test_rmse_mean[cv$best_iteration]

      # update best if improved
      if (rmse < best_rmse) {
        best_rmse    <- rmse
        best_params  <- params
        best_nrounds <- cv$best_iteration
      }
    }

    # close progress bar
    close(pb)

    # return best results
    list(
      params  = best_params,
      nrounds = best_nrounds,
      rmse    = best_rmse,
      grid_n  = nrow(grid)
    )
  }

  # message so it is clear why the simulation progress bar has not started yet
  cat("\nTuning the XGB model (one-time CV)...\n")
  cat(" - Profile:", tune_profile, "\n")
  cat(" - CV folds:", nfold, "\n")
  cat(" - Max grid size:", max_grid, "\n")
  cat(" - nrounds_max:", nrounds_max, " | early_stopping_rounds:", early_stopping_rounds, "\n\n")

  # tune X and Y separately
  # stack data for tuning``
  y_stack_x <- unlist(df[x_cols])
  C_stack_x <- C[rep(seq_len(nrow(C)), times = length(x_cols)), , drop = FALSE]

  # message for X tuning
  cat(" - Tuning X residualiser...\n")

  # tune X
  tuned_x <- tune_one(y_stack_x, C_stack_x, n_person = nrow(C))

  # stack data for tuning
  y_stack_y <- unlist(df[y_cols])
  C_stack_y <- C[rep(seq_len(nrow(C)), times = length(y_cols)), , drop = FALSE]

  # message for Y tuning
  cat("\n - Tuning Y residualiser...\n")

  # tune Y
  tuned_y <- tune_one(y_stack_y, C_stack_y, n_person = nrow(C))
  
  # newline so the next progress bar starts on a clean line
  cat("\n")

  # return results
  out <- list(
    tune_profile = tune_profile,

    params_X  = tuned_x$params,
    nrounds_X = tuned_x$nrounds,
    rmse_X    = tuned_x$rmse,
    grid_n_X  = tuned_x$grid_n,

    params_Y  = tuned_y$params,
    nrounds_Y = tuned_y$nrounds,
    rmse_Y    = tuned_y$rmse,
    grid_n_Y  = tuned_y$grid_n,

    nfold = nfold,
    nrounds_max = nrounds_max,
    early_stopping_rounds = early_stopping_rounds,
    max_grid = max_grid
  )

  # make tuned settings available to residualise_panel_xgb() in this R session
  XGB_TUNED <<- out

  return(out)
}
