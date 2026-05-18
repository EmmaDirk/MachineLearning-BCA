# Functions to safely fit models using lavaan, capturing errors without stopping execution.
# The BCA fitting functions call the residualisation helper function before fitting.
#
# The function fit the following models:
# - 1. CLPM without confounder adjustment
# - 2. CLPM with direct confounder adjustment
# - 3. RI-CLPM with indirect confounder adjustment via random intercepts
# - 4. DPM
# - 5. RI-CLPM with freed latent factor loadings
# - 6. DPM with freed latent factor loadings
# - 7. BCA CLPM with residualised X and Y
# - 8. BCA RI-CLPM with residualised X and Y
# - 9. BCA DPM with residualised X and Y
# ------------------------------------------------------------------------------------------------------------

# 1. CLPM without confounder adjustment
# ------------------------------------------------------------
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

# 2. CLPM with direct confounder adjustment
# -------------------------------------------------------------
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

# 3. RI-CLPM with indirect confounder adjustment via random intercepts
# -------------------------------------------------------------
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

# 4. DPM
# --------------------------------------------------------------
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

# 5. RI-CLPM with freed latent factor loadings
# --------------------------------------------------------------
safe_fit_riclpm_free_loadings <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "ML",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# 6. DPM with freed latent factor loadings
# --------------------------------------------------------------
safe_fit_dpm_free_loadings <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # try to fit
  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "ML",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# 7. CLPM with residualisation of X and Y (BCA CLPM)
# --------------------------------------------------------------
safe_fit_clpm_resid <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # first residualise the data
  df_resid <- tryCatch(

    # residualise the data using the helper function
    residualise_panel_linearC(data),

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

# 8. RI-CLPM with residualisation of X and Y (BCA RI-CLPM)
# --------------------------------------------------------------
safe_fit_riclpm_resid <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # first residualise the data
  df_resid <- tryCatch(
    residualise_panel_linearC(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # if residualisation failed, return NULL fit and the error message
  if (is.null(df_resid)) {
    return(list(fit = NULL, err = err))
  }

  # try to fit the RI-CLPM on the residualised data
  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_resid,
      estimator = "ML",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

# 9. DPM with residualisation of X and Y (BCA DPM)
# --------------------------------------------------------------
safe_fit_dpm_resid <- function(model_string, data) {

  # initialize error message
  err <- NA_character_

  # first residualise the data
  df_resid <- tryCatch(
    residualise_panel_linearC(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # if residualisation failed, return NULL fit and the error message
  if (is.null(df_resid)) {
    return(list(fit = NULL, err = err))
  }

  # try to fit the DPM on the residualised data
  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_resid,
      estimator = "ML",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  # return fit and error
  list(fit = fit, err = err)
}

