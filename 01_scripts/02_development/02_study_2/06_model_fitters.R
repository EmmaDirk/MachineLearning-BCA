# Functions to safely fit models using lavaan, capturing errors without stopping execution.
# Two BCA-style fitting functions first residualise X and Y on the confounders, then fit a CLPM to the residualised data:
# - one uses linear residualisation
# - one uses XGBoost (XGB) residualisation
#
# The functions fit the following models:
# - CLPM without confounder adjustment
# - CLPM with direct confounder adjustment (delta coefficients)
# - RI-CLPM with indirect confounder adjustment via random intercepts
# - DPM
# - BCA CLPM with linearly residualised X and Y
# - BCA CLPM with XGB residualised X and Y
#
# Naming logic:
# - autoregressive effects are beta, cross-lagged effects are gamma (these live in the lag matrix A elsewhere)
# - confounder effects are delta (these live in the delta matrices D_t elsewhere)
# - in the simulation output, these model-fitting functions are referenced by method labels like:
#   CLPM, RI-CLPM, DPM, CLPM_Adj, CLPM_LBCA (linear BCA), and CLPM_XGB (XGB BCA)
# ------------------------------------------------------------------------------------------------------------

safe_fit_clpm <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the data
      data      = as.data.frame(data),

      # use full information maximum likelihood
      estimator = "ML",
      
      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# same as above but for RI-CLPM
safe_fit_riclpm <- function(model_string, data) {
  
  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the data
      data      = as.data.frame(data),
      estimator = "ML",

      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# same as above but for DPM
safe_fit_dpm <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the data
      data      = as.data.frame(data),

      # use full information maximum likelihood
      estimator = "ML",

      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# same as above but for CLPM with confounders (delta coefficients)
safe_fit_clpm_C <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the data
      data      = as.data.frame(data),

      # use full information maximum likelihood
      estimator = "ML",

      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# same as above but for CLPM with linear BCA residualisation
safe_fit_clpm_resid <- function(model_string, data, k) {

  # initialize error message
  err <- NA_character_

  # first residualise the data
  df_resid <- tryCatch(

    # residualise the data using the helper function
    residualise_panel_linearC(data, k),

    # capture error message if residualisation fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # if residualisation failed, return NULL fit and the error message
  if (is.null(df_resid)) {
    return(list(fit = NULL, err = err))
  }

  # try to fit the CLPM on the residualised data
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the residualised data
      data      = df_resid,

      # use full information maximum likelihood
      estimator = "ML",

      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# same as above but for CLPM with XGB BCA residualisation
safe_fit_clpm_xgb <- function(model_string, data, k) {

  # initialize error message
  err <- NA_character_

  # first residualise the data using XGBoost
  df_resid <- tryCatch(

    # residualise the data using the XGBoost helper function
    residualise_panel_xgb(data, k),

    # capture error message if residualisation fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # if residualisation failed, return NULL fit and the error message
  if (is.null(df_resid)) {
    return(list(fit = NULL, err = err))
  }

  # try to fit the CLPM on the residualised data
  fit <- tryCatch(

    # use lavaan
    lavaan::lavaan(

      # the model string produced by the model builder
      model_string,

      # the residualised data
      data      = df_resid,

      # use full information maximum likelihood
      estimator = "ML",

      # turn off warnings
      warn      = FALSE
    ),

    # capture error message if fitting fails
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}
