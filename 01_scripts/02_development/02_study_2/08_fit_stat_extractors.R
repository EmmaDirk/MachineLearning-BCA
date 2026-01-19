# since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive, cross-lagged, and residual correlation (rho) are located
# this file contains functions to extract those parameters from the fitted model objects
# ------------------------------------------------------------------------------------------------

# lagged parameters extractor
# returns: est, p-value, CI lower, CI upper for each lagged parameter
extract_lagged_parameters <- function(
    fit,                                                      # lavaan model object
    T,                                                        # number of time points
    model_type = c("clpm", "riclpm", "dpm"),                  # model type
    ci_level = 0.95                                          # CI level
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
    ci_level = 0.95                                                # (NEW) CI level for rho
){

  # match model type
  model_type <- match.arg(model_type)

  # (NEW) return a stats table for rho: est, p, CI lower/upper
  empty_rho <- data.frame(
    est      = rep(NA_real_, T),
    p        = rep(NA_real_, T),
    ci.lower = rep(NA_real_, T),
    ci.upper = rep(NA_real_, T)
  )

  # if the model fit failed, return NAs
  if (is.null(fit)) return(empty_rho)

  # (NEW) use standardizedSolution so we directly get *correlations* with SE/p/CI
  ss <- tryCatch(
    lavaan::standardizedSolution(fit, se = TRUE, ci = TRUE, level = ci_level),
    error = function(e) NULL
  )

  # if extraction failed, return NAs
  if (is.null(ss)) return(empty_rho)

  # determine variable names based on model type, default is x, otherwise is wx
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # prepare container
  rho <- empty_rho

  # extract rho
  for (t in 1:T) {

    # the left hand side of the correlation equation
    lhs_xy <- paste0(xvar, t)

    # the right hand side of the correlation equation
    lhs_yx <- paste0(yvar, t)

    # (NEW) find the standardized covariance entry (this is a correlation in standardizedSolution)
    ix <- which(ss$op == "~~" & ss$lhs == lhs_xy & ss$rhs == lhs_yx)

    # if not found, try the other direction
    if (length(ix) == 0) {
      ix <- which(ss$op == "~~" & ss$lhs == lhs_yx & ss$rhs == lhs_xy)
    }

    # if not found, return NA
    if (length(ix) == 0) {
      rho$est[t]      <- NA_real_
      rho$p[t]        <- NA_real_
      rho$ci.lower[t] <- NA_real_
      rho$ci.upper[t] <- NA_real_
      next
    }

    i1 <- ix[1]

    # (NEW) standardized estimate column is typically est.std
    est <- if ("est.std" %in% names(ss)) ss$est.std[i1] else NA_real_

    # (NEW) pvalue + CI, when se=TRUE, ci=TRUE
    pval <- if ("pvalue"   %in% names(ss)) ss$pvalue[i1]   else NA_real_
    cil  <- if ("ci.lower" %in% names(ss)) ss$ci.lower[i1] else NA_real_
    ciu  <- if ("ci.upper" %in% names(ss)) ss$ci.upper[i1] else NA_real_

    rho$est[t]      <- est
    rho$p[t]        <- pval
    rho$ci.lower[t] <- cil
    rho$ci.upper[t] <- ciu
  }

  # return the residual correlations vector
  rho
}