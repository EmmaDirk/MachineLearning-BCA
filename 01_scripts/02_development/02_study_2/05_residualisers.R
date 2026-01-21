# For BCA (Baseline Covariate Adjustment) in panel data, we need to residualise each observed variable at each time point.
# Concretely, we fit the models: x_t ~ f(U) + e, and y_t ~ f(U) + e, where U are the baseline confounders with linear effects. 
# We then replace the observed x_t with e (the residuals), and similarly for y_t. 
#
# We use two functions for F(U). 
# - residualise_panel_linearC() uses a linear regression model (only using the linear confounders)
# - residualise_panel_xgb() uses an XGBoost model (using only the linear confounders, but can capture non-linear relationships)
# --------------------------------------------------------------------------------

# linear residualiser
residualise_panel_linearC <- function(df,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  
  # convert to data frame
  df <- as.data.frame(df)
  
  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value=TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value=TRUE)
  c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value=TRUE)

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # convert confounders to matrix
  C <- as.matrix(df[c_cols])

  # for each x and y, residualise against confounders
  for (x in x_cols)

    # with the linear model: x_t ~ confounders, and replace the column with the residuals
    df[[x]] <- resid(lm(df[[x]] ~ C))

  # same for y
  for (y in y_cols)
    df[[y]] <- resid(lm(df[[y]] ~ C))

  # return the residualised data frame
  df
}

# XGB residualiser
residualise_panel_xgb <- function(df,
                                  k,
                                  x_prefix = "x",
                                  y_prefix = "y",
                                  c_prefix = "c") {

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

  # confounder matrix (standardised for XGB stability)
  C <- scale(as.matrix(df[c_cols]))

  # pull tuned settings
  if (!exists("XGB_TUNED", inherits = TRUE) || is.null(XGB_TUNED))
    stop("XGB_TUNED not found. Tune once before calling residualise_panel_xgb().")

  # expect separate tuned settings for X and Y
  # i.e. we use double xgb tuning: once for X's, once for Y's. If one model is correct, we are already good. 
  if (is.null(XGB_TUNED$params_X) || is.null(XGB_TUNED$nrounds_X) ||
      is.null(XGB_TUNED$params_Y) || is.null(XGB_TUNED$nrounds_Y)) {
    stop("XGB_TUNED must contain params_X/nrounds_X and params_Y/nrounds_Y. Re-run tune_xgb_once().")
  }

  # pull tuned settings
  # for X's
  nrounds_X <- XGB_TUNED$nrounds_X
  params_X  <- XGB_TUNED$params_X

  # for Y's
  nrounds_Y <- XGB_TUNED$nrounds_Y
  params_Y  <- XGB_TUNED$params_Y

  # helper: fit xgboost and return residuals
  xgb_resid <- function(y, params_tuned, nrounds_tuned) {

    # create DMatrix
    dtrain <- xgboost::xgb.DMatrix(data = C, label = y)

    # fit
    fit <- xgboost::xgb.train(

      # data
      data = dtrain,

      # number of rounds and parameters
      nrounds = nrounds_tuned,

      # fixed + tuned parameters
      params = c(
        list(

          # objective and evaluation
          objective = "reg:squarederror",
          eval_metric = "rmse",

          # tree settings
          booster   = "gbtree",
          tree_method = "hist",

          # avoid oversubscribing CPU when the main simulation is parallel
          nthread = 1
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
  for (x in x_cols)
    df[[x]] <- xgb_resid(df[[x]], params_X, nrounds_X)

  # residualise y's
  for (y in y_cols)
    df[[y]] <- xgb_resid(df[[y]], params_Y, nrounds_Y)

  # return residualised data frame
  df
}