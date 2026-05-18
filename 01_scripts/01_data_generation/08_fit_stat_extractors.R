# =================================================================================================
# 
# These functions extract exactly the quantities needed for the final simulation output:
# - lagged estimates from the fitted SEM
# - model-based standard errors from the fitted SEM
# - the BIC
# - the true lagged parameters implied by Phi
#
# Everything is expanded to length T, with occasion 1 set to NA,
# because no lagged effect exists at the first occasion.
# =================================================================================================

# extract lagged estimates and their model-based standard errors
extract_lagged_estimates <- function(
    fit,
    T,
    model_type = c("clpm", "riclpm", "dpm")
) {

  # match model type
  model_type <- match.arg(model_type)

  # prepare empty vectors of length T; occasion 1 has no lagged effects
  out <- list(
    ARX    = rep(NA_real_, T),
    se_ARX = rep(NA_real_, T),
    ARY    = rep(NA_real_, T),
    se_ARY = rep(NA_real_, T),
    CXY    = rep(NA_real_, T),
    se_CXY = rep(NA_real_, T),
    CYX    = rep(NA_real_, T),
    se_CYX = rep(NA_real_, T)
  )

  # failed fit returns all NAs
  if (is.null(fit)) {
    return(out)
  }

  # extract parameter table
  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)
  if (is.null(pe)) {
    return(out)
  }

  # RI-CLPM uses within-person latent variables; the others use observed variables
  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  # helper to grab estimate and SE for one regression path
  grab <- function(lhs, rhs) {
    ix <- which(pe$op == "~" & pe$lhs == lhs & pe$rhs == rhs)
    if (length(ix) == 0) return(c(NA_real_, NA_real_))
    i <- ix[1]
    c(pe$est[i], if ("se" %in% names(pe)) pe$se[i] else NA_real_)
  }

  # fill all lagged effects from occasion 2 onward
  for (t in 2:T) {

    arx <- grab(paste0(xvar, t), paste0(xvar, t - 1))
    ary <- grab(paste0(yvar, t), paste0(yvar, t - 1))
    cxy <- grab(paste0(yvar, t), paste0(xvar, t - 1))
    cyx <- grab(paste0(xvar, t), paste0(yvar, t - 1))

    out$ARX[t]    <- arx[1]
    out$se_ARX[t] <- arx[2]
    out$ARY[t]    <- ary[1]
    out$se_ARY[t] <- ary[2]
    out$CXY[t]    <- cxy[1]
    out$se_CXY[t] <- cxy[2]
    out$CYX[t]    <- cyx[1]
    out$se_CYX[t] <- cyx[2]
  }

  out
}


# build the true lagged parameter vectors implied by Phi
extract_true_lagged_parameters <- function(Phi, T) {

  # basic check
  if (!is.matrix(Phi) || !all(dim(Phi) == c(2, 2))) {
    stop("Phi must be a 2 x 2 matrix.")
  }

  # occasion 1 has no lagged effects
  beta_x   <- rep(NA_real_, T)
  beta_y   <- rep(NA_real_, T)
  gamma_xy <- rep(NA_real_, T)
  gamma_yx <- rep(NA_real_, T)

  # constant lag matrix across occasions 2..T
  if (T >= 2) {
    beta_x[2:T]   <- Phi[1, 1]
    beta_y[2:T]   <- Phi[2, 2]
    gamma_xy[2:T] <- Phi[2, 1]
    gamma_yx[2:T] <- Phi[1, 2]
  }

  list(
    beta_x = beta_x,
    beta_y = beta_y,
    gamma_xy = gamma_xy,
    gamma_yx = gamma_yx
  )
}


# extract the BIC once and repeat it later across occasions
extract_bic <- function(fit) {

  if (is.null(fit)) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(lavaan::fitMeasures(fit, "bic")),
    error = function(e) NA_real_
  )
}
