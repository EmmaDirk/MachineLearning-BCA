# These helpers handle the parts of the workflow that are not the main SEM fit itself:
# - classifying runs as success / non-converged / improper
# - bootstrap-based standard errors for two-stage procedures
#
# Flag coding:
# - 0 = successful and proper run
# - 1 = failed or non-converged run
# - 2 = converged but improper solution
# -------------------------------------------------------------------------------------------------

# classify a fitted lavaan object into the requested flag coding
classify_fit_flag <- function(fit) {

  # completely failed fit counts as non-convergence
  if (is.null(fit)) {
    return(1L)
  }

  # check convergence first
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) {
    return(1L)
  }

  # lavaan already provides a broad admissibility check
  post_ok <- tryCatch(lavaan::lavInspect(fit, "post.check"), error = function(e) TRUE)
  if (!isTRUE(post_ok)) {
    return(2L)
  }

  # extra guard: catch negative variances if they appear in the parameter table
  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)
  if (!is.null(pe)) {
    var_rows <- pe$op == "~~" & pe$lhs == pe$rhs
    if (any(pe$est[var_rows] < 0, na.rm = TRUE)) {
      return(2L)
    }
  }

  # otherwise it is a successful run
  0L
}


# determine whether the chosen pipeline needs bootstrap-based standard errors
uses_bootstrap_se <- function(residualizer) {

  residualizer %in% c("linear", "xgb")
}


# bootstrap the whole analysis pipeline and return empirical SEs for the lagged effects
bootstrap_pipeline_se <- function(
    df,                                          # data frame
    T,                                           # number of time points
    k,                                           # number of confounders
    residualizer,                                # residualiser of choice
    sem_model,                                   # model of choice
    confounder_order,                            # order of confounders for residualisation and CLPM
    exclude,                                     # features to exclude from residualisation and CLPM
    free_loadings,                               # free loadings for DPM and RI-CLPM
    xgb_tuning,                                  # XGB hyperparameters
    residualizer_args,                           # residualiser-specific arguments
    B,                                           # number of bootstrap samples
    seed = NULL
) {

  # if B < 2, return empty SE vectors
  if (is.null(B) || B < 2) {
    return(list(
      ARX = rep(NA_real_, T),
      ARY = rep(NA_real_, T),
      CXY = rep(NA_real_, T),
      CYX = rep(NA_real_, T)
    ))
  }

  # optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  # ensure that each original row carries a stable id through the bootstrap
  if (!(".id_orig" %in% names(df))) {
    df$.id_orig <- seq_len(nrow(df))
  }

  # storage for bootstrap estimates; occasion 1 stays NA by design
  ARX <- matrix(NA_real_, nrow = B, ncol = T)
  ARY <- matrix(NA_real_, nrow = B, ncol = T)
  CXY <- matrix(NA_real_, nrow = B, ncol = T)
  CYX <- matrix(NA_real_, nrow = B, ncol = T)

  # bootstrap the full pipeline
  for (b in seq_len(B)) {

    # sample rows with replacement
    idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
    df_b <- df[idx, , drop = FALSE]

    # refit the entire pipeline on the resample
    fit_b <- fit_analysis_pipeline(
      df = df_b,
      T = T,
      k = k,
      residualizer = residualizer,
      sem_model = sem_model,
      confounder_order = confounder_order,
      exclude = exclude,
      free_loadings = free_loadings,
      xgb_tuning = xgb_tuning,
      residualizer_args = residualizer_args
    )

    # skip failed bootstrap fits
    if (is.null(fit_b$fit)) next

    # extract lagged estimates from the bootstrap fit
    lag_b <- extract_lagged_estimates(
      fit = fit_b$fit,
      T = T,
      model_type = sem_model
    )

    ARX[b, ] <- lag_b$ARX
    ARY[b, ] <- lag_b$ARY
    CXY[b, ] <- lag_b$CXY
    CYX[b, ] <- lag_b$CYX
  }

  # empirical standard errors across bootstrap replications
  list(
    ARX = apply(ARX, 2, stats::sd, na.rm = TRUE),
    ARY = apply(ARY, 2, stats::sd, na.rm = TRUE),
    CXY = apply(CXY, 2, stats::sd, na.rm = TRUE),
    CYX = apply(CYX, 2, stats::sd, na.rm = TRUE)
  )
}