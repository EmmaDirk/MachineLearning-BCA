# 07_bootstrap_helpers.R
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


# convert one integer flag into three proportions that sum to 1
flag_to_props <- function(flag) {

  out <- c(flag0 = 0, flag1 = 0, flag2 = 0)

  if (is.na(flag) || !(flag %in% c(0L, 1L, 2L, 0, 1, 2))) {
    return(as.list(out))
  }

  out[paste0("flag", as.integer(flag))] <- 1

  as.list(out)
}


# helper to turn a variance name into a readable diagnosis
describe_negative_variance_target <- function(name) {

  if (is.na(name) || !nzchar(name)) {
    return("negative variance")
  }

  # RI-CLPM random intercepts
  if (grepl("^ri[xy]$", name)) {
    return(sprintf("negative random intercept variance (%s)", name))
  }

  # DPM accumulating factors
  if (grepl("^F[XY]$", name)) {
    return(sprintf("negative accumulating factor variance (%s)", name))
  }

  # RI-CLPM within-person factors
  if (grepl("^w[xy][0-9]+$", name)) {
    return(sprintf("negative within-person latent variance (%s)", name))
  }

  # observed variables / residual variances
  if (grepl("^[xy][0-9]+$", name)) {
    return(sprintf("negative observed residual variance (%s)", name))
  }

  # generic latent or observed variance
  sprintf("negative variance (%s)", name)
}


# helper to detect whether a covariance matrix is not positive semidefinite
matrix_not_psd <- function(S, tol = 1e-10) {

  if (is.null(S)) {
    return(FALSE)
  }

  if (!is.matrix(S)) {
    S <- tryCatch(as.matrix(S), error = function(e) NULL)
  }

  if (is.null(S)) {
    return(FALSE)
  }

  if (nrow(S) == 0 || ncol(S) == 0) {
    return(FALSE)
  }

  ev <- tryCatch(eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)

  if (all(is.na(ev))) {
    return(FALSE)
  }

  any(ev < -tol, na.rm = TRUE)
}


# diagnose the main reason why a converged fit is improper
diagnose_improper_fit <- function(fit, tol = 1e-10) {

  # failed fit
  if (is.null(fit)) {
    return(NA_character_)
  }

  # non-converged fit
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) {
    return(NA_character_)
  }

  # if lavaan thinks the fit is admissible, do not attach a reason
  post_ok <- tryCatch(lavaan::lavInspect(fit, "post.check"), error = function(e) TRUE)

  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)

  # 1) first and most interpretable case: explicit negative variances
  if (!is.null(pe)) {

    var_rows <- pe$op == "~~" & pe$lhs == pe$rhs & !is.na(pe$est)

    if (any(var_rows)) {

      neg_var_rows <- which(var_rows & pe$est < -tol)

      if (length(neg_var_rows) > 0) {

        # choose the most negative variance as the main culprit
        i <- neg_var_rows[which.min(pe$est[neg_var_rows])]
        return(describe_negative_variance_target(pe$lhs[i]))
      }
    }
  }

  # 2) check the latent covariance matrix
  cov_lv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
  if (matrix_not_psd(cov_lv, tol = tol)) {
    return("latent covariance matrix not positive semidefinite")
  }

  # 3) check the observed residual covariance matrix
  theta <- tryCatch(lavaan::lavInspect(fit, "theta"), error = function(e) NULL)
  if (matrix_not_psd(theta, tol = tol)) {
    return("observed residual covariance matrix not positive semidefinite")
  }

  # 4) if lavaan flagged the solution as improper but we did not localize it more precisely
  if (!isTRUE(post_ok)) {
    return("post.check failed: unspecified improper solution")
  }

  # otherwise there is no improper-fit reason to report
  NA_character_
}


# convert fit status into a compact bootstrap issue label
diagnose_bootstrap_issue <- function(fit, tol = 1e-10) {

  flag <- classify_fit_flag(fit)

  if (flag == 0L) {
    return("proper")
  }

  if (flag == 1L) {
    return("nonconverged_or_failed")
  }

  reason <- diagnose_improper_fit(fit, tol = tol)

  if (is.na(reason) || !nzchar(reason)) {
    return("improper_unspecified")
  }

  reason
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
      CYX = rep(NA_real_, T),
      bootstrap_prop_success = NA_real_,
      flag0 = NA_real_,
      flag1 = NA_real_,
      flag2 = NA_real_,
      bootstrap_issue_vector = NA_character_
    ))
  }

  # if no bootstrap seed was supplied, generate one automatically
  # this keeps the function usable even when seed = NULL
  if (is.null(seed)) {
    max_seed <- max(1L, .Machine$integer.max - as.integer(B) - 1L)
    seed <- as.integer(sample.int(max_seed, size = 1))
  } else {
    seed <- as.integer(seed[1])
  }

  # set the bootstrap seed once
  set.seed(seed)

  # ensure that each original row carries a stable id through the bootstrap
  if (!(".id_orig" %in% names(df))) {
    df$.id_orig <- seq_len(nrow(df))
  }

  # storage for bootstrap estimates; occasion 1 stays NA by design
  ARX <- matrix(NA_real_, nrow = B, ncol = T)
  ARY <- matrix(NA_real_, nrow = B, ncol = T)
  CXY <- matrix(NA_real_, nrow = B, ncol = T)
  CYX <- matrix(NA_real_, nrow = B, ncol = T)

  # track the bootstrap flag on every resample
  boot_flag <- rep(NA_integer_, B)

  # track the bootstrap issue on every resample
  bootstrap_issue_vector <- rep(NA_character_, B)

  # bootstrap the full pipeline
  for (b in seq_len(B)) {

    # sample rows with replacement
    idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
    df_b <- df[idx, , drop = FALSE]

    # vary the stage-1 fold allocation seed across bootstrap draws
    # if we keep the same fold seed in every draw, the bootstrap misses an
    # important source of stage-1 uncertainty and the SEs can become too small
    residualizer_args_b <- residualizer_args

    if (is.null(residualizer_args_b)) {
      residualizer_args_b <- list()
    }

    residualizer_args_b$seed <- as.integer(seed + b)

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
      residualizer_args = residualizer_args_b
    )

    # classify every bootstrap fit, including outright failures
    boot_flag[b] <- classify_fit_flag(fit_b$fit)

    # store a readable bootstrap issue label
    bootstrap_issue_vector[b] <- diagnose_bootstrap_issue(fit_b$fit)

    # skip failed bootstrap fits for estimate extraction
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
    CYX = apply(CYX, 2, stats::sd, na.rm = TRUE),

    # keep the old success metric for compatibility
    bootstrap_prop_success = mean(boot_flag == 0L, na.rm = TRUE),

    # proportional flag outputs
    flag0 = mean(boot_flag == 0L, na.rm = TRUE),
    flag1 = mean(boot_flag == 1L, na.rm = TRUE),
    flag2 = mean(boot_flag == 2L, na.rm = TRUE),

    # full bootstrap issue trace
    bootstrap_issue_vector = bootstrap_issue_vector
  )
}