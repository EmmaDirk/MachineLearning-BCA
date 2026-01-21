# For BCA (Baseline Covariate Adjustment) in panel data, we need to residualise each observed variable at each time point.
# Concretely, we fit the models: x_t ~ f(U) + e, and y_t ~ f(U) + e, where U are the baseline confounders with linear effects. 
# We then replace the observed x_t with e (the residuals), and similarly for y_t. 
#
# We use two functions for F(U). 
# - residualise_panel_linearC() uses a linear regression model (only using the linear confounders)
# - residualise_panel_xgb() uses an XGBoost model (using only the linear confounders, but can capture non-linear relationships)
# --------------------------------------------------------------------------------

# new xgb tuning results are passed in as an argument
# new xgb fit speed can be controlled by a fit profile argument

# linear residualiser
residualise_panel_linearC <- function(df,
                                      k = NULL,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  
  # convert to data frame
  df <- as.data.frame(df)
  
  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value=TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value=TRUE)

  # choose confounders
  if (is.null(k)) {

    # use all confounder columns that match the prefix
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value=TRUE)

  } else {

    # use only linear confounders c1 to ck
    c_cols <- paste0(c_prefix, 1:k)
  }

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0)
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))

  # convert confounders to matrix
  C <- as.matrix(df[c_cols])

  # for each x and y, residualise against confounders
  for (x in x_cols) {

    # with the linear model: x_t ~ confounders, and replace the column with the residuals
    df[[x]] <- resid(lm(df[[x]] ~ C))
  }

  # same for y
  for (y in y_cols) {
    df[[y]] <- resid(lm(df[[y]] ~ C))
  }

  # return the residualised data frame
  df
}

# XGB residualiser
residualise_panel_xgb <- function(df,
                                  k,
                                  xgb_tuned,
                                  fit_profile = c("fast", "balanced", "thorough"),
                                  fit_overrides = NULL,
                                  x_prefix = "x",
                                  y_prefix = "y",
                                  c_prefix = "c") {

  # match fit profile argument
  fit_profile <- match.arg(fit_profile)

  # convert overrides to list
  if (is.null(fit_overrides)) fit_overrides <- list()
  if (!is.list(fit_overrides)) stop("fit_overrides must be a list or NULL")

  # convert to data frame
  df <- as.data.frame(df)

  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # only linear confounders are observed
  c_cols <- paste0(c_prefix, 1:k)

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if requested confounders are missing
  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0)
    stop("Requested confounder columns not found: ", paste(missing_c, collapse = ", "))

  # confounder matrix (standardised for XGB stability)
  C <- scale(as.matrix(df[c_cols]))

  # pull tuned settings
  if (missing(xgb_tuned) || is.null(xgb_tuned))
    stop("xgb_tuned not provided. Tune once before calling residualise_panel_xgb().")

  # expect separate tuned settings for X and Y
  # i.e. we use double xgb tuning: once for X's, once for Y's. If one model is correct, we are already good. 
  if (is.null(xgb_tuned$params_X) || is.null(xgb_tuned$nrounds_X) ||
      is.null(xgb_tuned$params_Y) || is.null(xgb_tuned$nrounds_Y)) {
    stop("xgb_tuned must contain params_X nrounds_X params_Y and nrounds_Y")
  }

  # pull tuned settings
  # for X's
  nrounds_X <- xgb_tuned$nrounds_X
  params_X  <- xgb_tuned$params_X

  # for Y's
  nrounds_Y <- xgb_tuned$nrounds_Y
  params_Y  <- xgb_tuned$params_Y

  # defaults for fit speed
  fit_defaults <- switch(
    fit_profile,
    fast = list(
      nthread = 1,
      tree_method = "hist",
      max_bin = 128,
      nrounds_factor = 0.5
    ),
    balanced = list(
      nthread = 1,
      tree_method = "hist",
      max_bin = 256,
      nrounds_factor = 1.0
    ),
    thorough = list(
      nthread = 1,
      tree_method = "hist",
      max_bin = 512,
      nrounds_factor = 1.0
    )
  )

  # apply overrides
  if (!is.null(fit_overrides$nthread)) fit_defaults$nthread <- fit_overrides$nthread
  if (!is.null(fit_overrides$tree_method)) fit_defaults$tree_method <- fit_overrides$tree_method
  if (!is.null(fit_overrides$max_bin)) fit_defaults$max_bin <- fit_overrides$max_bin
  if (!is.null(fit_overrides$nrounds_factor)) fit_defaults$nrounds_factor <- fit_overrides$nrounds_factor

  # helper to compute the number of rounds to use
  nrounds_use <- function(nrounds_tuned) {

    # allow explicit override
    if (!is.null(fit_overrides$nrounds)) {
      return(max(1, as.integer(fit_overrides$nrounds)))
    }

    # otherwise scale tuned rounds
    max(1, as.integer(floor(nrounds_tuned * fit_defaults$nrounds_factor)))
  }

  # helper: fit xgboost and return residuals
  xgb_resid <- function(y, params_tuned, nrounds_tuned) {

    # create DMatrix
    dtrain <- xgboost::xgb.DMatrix(data = C, label = y)

    # number of rounds to use
    nrounds_now <- nrounds_use(nrounds_tuned)

    # fit
    fit <- xgboost::xgb.train(

      # data
      data = dtrain,

      # number of rounds and parameters
      nrounds = nrounds_now,

      # fixed plus tuned parameters
      params = c(
        list(

          # objective and evaluation
          objective = "reg:squarederror",
          eval_metric = "rmse",

          # tree settings
          booster   = "gbtree",
          tree_method = fit_defaults$tree_method,
          max_bin = fit_defaults$max_bin,

          # avoid oversubscribing CPU when the main simulation is parallel
          nthread = fit_defaults$nthread
        ),

        # tuned parameters
        params_tuned
      ),

      # silent
      verbose = 0
    )

    # predicted values
    yhat <- predict(fit, dtrain)

    # return residuals
    y - yhat
  }

  # residualise x's
  for (x in x_cols) {
    df[[x]] <- xgb_resid(df[[x]], params_X, nrounds_X)
  }

  # residualise y's
  for (y in y_cols) {
    df[[y]] <- xgb_resid(df[[y]], params_Y, nrounds_Y)
  }

  # return residualised data frame
  df
}
