# Since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive (betas), cross-lagged (gammas), and residual correlation (rho) are located
# 
# Currently, the functions are equiped to extract parameters from:
# - CLPM
# - RI-CLPM
# - DPM
# ------------------------------------------------------------------------------------------------

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
    model_type = c("clpm","riclpm","dpm"),                         # model type
    ci_level = 0.95                                                # CI level
){

  # match model type
  model_type <- match.arg(model_type)

  # prepare empty output
  out <- data.frame(
    est = rep(NA_real_, T),
    p = rep(NA_real_, T),
    ci.lower = rep(NA_real_, T),
    ci.upper = rep(NA_real_, T)
  )

  # if the model fit failed, return NAs
  if (is.null(fit)) return(out)

  # try to extract parameter estimates with ci
  pe <- tryCatch(
    lavaan::parameterEstimates(fit, ci = TRUE, level = ci_level),
    error = function(e) NULL
  )

  # if extraction failed, return NAs
  if (is.null(pe)) return(out)

  # determine variable names based on model type, default is x, otherwise is wx
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # extract rho
  for (t in 1:T) {

    lhs_xy <- paste0(xvar, t)
    lhs_yx <- paste0(yvar, t)

    # covariance row
    ix <- which(pe$op == "~~" & pe$lhs == lhs_xy & pe$rhs == lhs_yx)
    if (length(ix) == 0) {
      ix <- which(pe$op == "~~" & pe$lhs == lhs_yx & pe$rhs == lhs_xy)
    }

    # variance rows
    vx_idx <- which(pe$op == "~~" & pe$lhs == lhs_xy & pe$rhs == lhs_xy)
    vy_idx <- which(pe$op == "~~" & pe$lhs == lhs_yx & pe$rhs == lhs_yx)

    if (length(ix) == 0 || length(vx_idx) == 0 || length(vy_idx) == 0) {
      next
    }

    cov_xy <- pe$est[ix[1]]
    vx <- pe$est[vx_idx[1]]
    vy <- pe$est[vy_idx[1]]

    if (is.na(vx) || is.na(vy) || vx <= 0 || vy <= 0) {
      next
    }

    out$est[t] <- cov_xy / sqrt(vx * vy)

    # p values and ci for rho are not returned directly by lavaan
    # keep these as NA to avoid misleading inference
  }

  out
}
