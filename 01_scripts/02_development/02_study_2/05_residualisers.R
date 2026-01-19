# function that residualises replaces x and y columns in a data frame
# by their residuals after regressing out confounders c1...ck
# --------------------------------------------------------------------------------

# this function is doing the exact same as a linearly adjusted CLPM, but instead of including the confounders in the model,
# we decouple the confounder adjustment from the model fitting by residualising all X and Y variables against the confounders
# this is called Baseline Covariate Adjustment (BCA)
residualise_panel_linearC <- function(df,
                                      k,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  
  # convert to data frame
  df <- as.data.frame(df)
  
  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # IMPORTANT: only use linear confounders (c1 ... ck)
  c_cols <- paste0(c_prefix, 1:k)

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # stop if some confounders are missing
  if (any(!c_cols %in% names(df)))
    stop("Not all linear confounder columns found (c1..ck).")

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

# now we want to also add a model that can deal with the non-linear relationships
# between the confounders and the outcome variables. This will be an Extreme Gradient Boosting Xgb model
# note that the model still only 'sees' the linear confounders, but since those are deterministically related 
# to the non-linear confounders, the Xgb model can in theory learn these non-linear relationships
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

  # (CHANGED) expect separate tuned settings for X and Y
  if (is.null(XGB_TUNED$params_X) || is.null(XGB_TUNED$nrounds_X) ||
      is.null(XGB_TUNED$params_Y) || is.null(XGB_TUNED$nrounds_Y)) {
    stop("XGB_TUNED must contain params_X/nrounds_X and params_Y/nrounds_Y. Re-run tune_xgb_once().")
  }

  nrounds_X <- XGB_TUNED$nrounds_X
  params_X  <- XGB_TUNED$params_X

  nrounds_Y <- XGB_TUNED$nrounds_Y
  params_Y  <- XGB_TUNED$params_Y

  # helper: fit xgboost and return residuals
  xgb_resid <- function(y, params_tuned, nrounds_tuned) {

    # create DMatrix
    dtrain <- xgboost::xgb.DMatrix(data = C, label = y)

    # fit
    fit <- xgboost::xgb.train(
      data = dtrain,
      nrounds = nrounds_tuned,
      params = c(
        list(
          objective = "reg:squarederror",
          eval_metric = "rmse",
          booster   = "gbtree",
          tree_method = "hist",

          # avoid oversubscribing CPU when the main simulation is parallel
          nthread = 1
        ),
        params_tuned
      ),
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