# Since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive (betas), cross-lagged (gammas), and residual correlation (rho) are located
# 
# Currently, the functions are equiped to extract parameters from (multiple versions of):
# - CLPM
# - RI-CLPM
# - DPM
# 
# this script additionally checks if models converge and if the solutions are proper
# i.e.:
# - no negative variances in any matrix (e.g. psi, theta, sigma)
# - all parameters are finite (no Inf or NaN)
# - all matrices are positive semidefinite
# ------------------------------------------------------------------------------------------------

# helper function to check for positive semidefiniteness
is_psd <- function(M, tol = 1e-8) {

  # return TRUE for NULL or non-matrix input (nothing to check)
  if (is.null(M) || !is.matrix(M)) return(TRUE)

  # 0x0 or 1x1 matrices are always psd
  if (nrow(M) <= 1) return(TRUE)

  # symmetrize defensively
  Ms <- (M + t(M)) / 2

  # eigenvalues
  ev <- tryCatch(
    eigen(Ms, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  # if eigen failed, we can't confirm psd, so return FALSE
  if (all(is.na(ev))) return(FALSE)

  # psd if smallest eigenvalue is not too negative
  min(ev) >= -tol
}

# helper function to check if the solution is proper
# returns a list with:
# - converged: did lavaan converge?
# - proper: is the solution proper?
# - reasons: character vector with reasons if not proper
check_convergence_and_properness <- function(fit, tol = 1e-8) {

  # if no fit, return NAs
  if (is.null(fit)) {
    return(list(
      converged = FALSE,
      proper    = NA,
      reasons   = NA_character_
    ))
  }

  # check convergence
  converged <- isTRUE(tryCatch(lavaan::inspect(fit, "converged"), error = function(e) FALSE))

  # if not converged, do not classify properness
  if (!converged) {
    return(list(
      converged = FALSE,
      proper    = NA,
      reasons   = "Model did not converge"
    ))
  }

  # collect reasons for improperness
  reasons <- character(0)

  # try to extract parameter estimates (also used for negative variances / finiteness)
  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)

  # check parameter finiteness
  if (is.null(pe)) {

    reasons <- c(reasons, "Could not extract parameter estimates")

  } else {

    # check for non-finite estimates
    if (any(!is.finite(pe$est))) {
      reasons <- c(reasons, "Non-finite parameter estimate(s) (Inf/NaN/NA)")
    }

    # check for non-finite standard errors if available
    if ("se" %in% names(pe) && any(!is.finite(pe$se))) {
      reasons <- c(reasons, "Non-finite standard error(s) (Inf/NaN/NA)")
    }

    # check for negative variances (free parameters only if possible)
    var_rows <- pe$op == "~~" & pe$lhs == pe$rhs

    # if free column exists, only consider freely estimated variances
    if ("free" %in% names(pe)) {
      var_rows <- var_rows & pe$free > 0
    }

    neg_vars <- pe[var_rows & is.finite(pe$est) & pe$est < -tol, , drop = FALSE]
    if (nrow(neg_vars) > 0) {
      reasons <- c(
        reasons,
        paste0("Negative variance(s): ", paste(unique(neg_vars$lhs), collapse = ", "))
      )
    }
  }

  # check key matrices for psd (theta, psi, sigma, cov.lv)
  theta  <- tryCatch(lavaan::lavInspect(fit, "theta"),  error = function(e) NULL)
  psi    <- tryCatch(lavaan::lavInspect(fit, "psi"),    error = function(e) NULL)
  sigma  <- tryCatch(lavaan::lavInspect(fit, "sigma"),  error = function(e) NULL)
  cov_lv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)

  if (!is_psd(theta,  tol = tol)) reasons <- c(reasons, "Theta (residual cov matrix) is not psd")
  if (!is_psd(psi,    tol = tol)) reasons <- c(reasons, "Psi (latent cov matrix) is not psd")
  if (!is_psd(sigma,  tol = tol)) reasons <- c(reasons, "Sigma (implied cov matrix) is not psd")
  if (!is_psd(cov_lv, tol = tol)) reasons <- c(reasons, "cov.lv (latent cov matrix) is not psd")

  # singularity / vcov problems (information matrix issues)
  V <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)

  if (is.null(V)) {
    reasons <- c(reasons, "vcov(fit) could not be computed (possible singular information matrix)")
  } else if (any(!is.finite(V))) {
    reasons <- c(reasons, "vcov(fit) contains non-finite values (possible singular information matrix)")
  }

  # proper if no reasons
  proper <- length(reasons) == 0

  list(
    converged = TRUE,
    proper    = proper,
    reasons   = if (!proper) reasons else character(0)
  )
}

# lagged parameters extractor
# now returns: est, p-value, CI lower, CI upper for each lagged parameter
extract_lagged_parameters <- function(
    fit,                                                      # lavaan model object
    T,                                                        # number of time points
    model_type = c("clpm", "riclpm", "dpm"),                  # model type
    ci_level = 0.95                                           # CI level
){

  # match model type
  model_type <- match.arg(model_type)

  # helper: empty stats data frame
  empty_stats <- data.frame(
    est      = rep(NA_real_, T-1),
    p        = rep(NA_real_, T-1),
    ci.lower = rep(NA_real_, T-1),
    ci.upper = rep(NA_real_, T-1)
  )

  # if the model fit failed, return NAs
  if (is.null(fit)) {
    return(list(
      ar_x = empty_stats,    # autoregressive X_t ← X_{t-1}
      ar_y = empty_stats,    # autoregressive Y_t ← Y_{t-1}
      xy   = empty_stats,    # cross-lag Y_t ← X_{t-1}   (X → Y)
      yx   = empty_stats     # cross-lag X_t ← Y_{t-1}   (Y → X)
    ))
  }

  # try to extract parameter table (with CI)
  pe <- tryCatch(
    lavaan::parameterEstimates(fit, ci = TRUE, level = ci_level),
    error = function(e) NULL
  )

  # if extraction failed, return NAs
  if (is.null(pe)) {
    return(list(
      ar_x = empty_stats,
      ar_y = empty_stats,
      xy   = empty_stats,
      yx   = empty_stats
    ))
  }

  # RI-CLPM uses latent within-person variables wx, wy
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # helper to grab a single parameter (regression paths only)
  grab_stats <- function(lhs, rhs) {

    ix <- which(pe$op == "~" & pe$lhs == lhs & pe$rhs == rhs)
    if (length(ix) == 0) return(c(NA_real_, NA_real_, NA_real_, NA_real_))

    i1 <- ix[1]

    est <- pe$est[i1]

    # p-value column is typically "pvalue"
    pval <- if ("pvalue" %in% names(pe)) pe$pvalue[i1] else NA_real_

    # CI columns are "ci.lower" and "ci.upper" when ci=TRUE
    cil <- if ("ci.lower" %in% names(pe)) pe$ci.lower[i1] else NA_real_
    ciu <- if ("ci.upper" %in% names(pe)) pe$ci.upper[i1] else NA_real_

    c(est, pval, cil, ciu)
  }

  # containers
  ar_x <- empty_stats
  ar_y <- empty_stats
  xy   <- empty_stats
  yx   <- empty_stats

  # extract all lagged parameters
  for (t in 2:T) {

    # autoregressive
    sx <- grab_stats(paste0(xvar, t), paste0(xvar, t-1))
    sy <- grab_stats(paste0(yvar, t), paste0(yvar, t-1))

    ar_x[t-1, ] <- sx
    ar_y[t-1, ] <- sy

    # cross-lag (X → Y)
    sxy <- grab_stats(paste0(yvar, t), paste0(xvar, t-1))
    xy[t-1, ] <- sxy

    # cross-lag (Y → X)
    syx <- grab_stats(paste0(xvar, t), paste0(yvar, t-1))
    yx[t-1, ] <- syx
  }

  list(
    ar_x = ar_x,
    ar_y = ar_y,
    xy   = xy,
    yx   = yx
  )
}

# residual correclations extractor
extract_rho_vec <- function(
    fit,                                                           # lavaan model object
    T,                                                             # number of time points
    model_type = c("clpm","riclpm","dpm")                          # model type
){

  # match model type
  model_type <- match.arg(model_type)

  # if the model fit failed, return NAs
  if (is.null(fit)) return(rep(NA_real_, T))

  # try to extract parameter estimates
  pe <- tryCatch(lavaan::parameterEstimates(fit), error=function(e) NULL)

  # if extraction failed, return NAs
  if (is.null(pe)) return(rep(NA_real_, T))

  # determine variable names based on model type, default is x, otherwise is wx
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # prepare container
  rho <- numeric(T)

  # extract rho
  for (t in 1:T) {

    # the left hand side of the correlation equation
    lhs_xy <- paste0(xvar, t)

    # the right hand side of the correlation equation
    lhs_yx <- paste0(yvar, t)

    # find the covariance estimate between x_t and y_t
    ix <- which(pe$op == "~~" & pe$lhs == lhs_xy & pe$rhs == lhs_yx)

    # if not found, try the other direction
    if (length(ix) == 0) {

      # find the covariance estimate between y_t and x_t
      ix <- which(pe$op == "~~" & pe$lhs == lhs_yx & pe$rhs == lhs_xy)
    }

    # if not found, return NA
    if (length(ix) == 0) {
      rho[t] <- NA_real_
      next
    }

    # covariance estimate
    cov_xy <- pe$est[ix[1]]

    # find the variance estimates for x_t and y_t
    vx_idx <- which(pe$op == "~~" & pe$lhs == lhs_xy & pe$rhs == lhs_xy)
    vy_idx <- which(pe$op == "~~" & pe$lhs == lhs_yx & pe$rhs == lhs_yx)

    # if not found, return NA
    if (length(vx_idx) == 0 || length(vy_idx) == 0) {
      rho[t] <- NA_real_
      next
    }

    # variance estimates
    vx <- pe$est[vx_idx[1]]
    vy <- pe$est[vy_idx[1]]

    # compute rho: cov_xy / sqrt(vx * vy)
    if (is.na(vx) || is.na(vy) || vx <= 0 || vy <= 0) {
      rho[t] <- NA_real_
    } else {
      rho[t] <- cov_xy / sqrt(vx * vy)
    }
  }

  # return the residual correlations vector
  rho
}
