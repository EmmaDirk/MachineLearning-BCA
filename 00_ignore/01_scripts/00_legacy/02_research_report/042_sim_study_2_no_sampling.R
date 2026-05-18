############################################################
# Simulation Study 2 Functions
############################################################
# 
# The goal of this script is to provide evidence that the RI-CLPM, just like the DPM, are not equipped
# to handle time-varying effects of baseline confounders. In this script we do the following:
# 1) (CHANGED) Remove sampling of B coefficients. Instead, allow the user to pass in B-matrices directly.
# 2) (CHANGED) Remove internal B trajectory generators (constant/stepwise). Instead, use user-supplied B lists per scenario.
# 3) We simulate panel data under a CLPM data-generating process with time-varying B matrices.
# 4) Build the models dynamically for adaption to T: CLPM, RI-CLPM, DPM, CLPM with confounders, and linear BCA in conjunction with CLPM. 
# 5) (CHANGED) For all extracted lagged parameters (AR and cross-lags), we now also return p-values and confidence intervals.

library(mvtnorm)
library(lavaan)
library(tidyverse)   
library(parallel)
library(pbapply)
library(xgboost)

############################################################
##  1. User-supplied B-matrices helpers (replaces sampling)
############################################################

# helper: compute number of interaction terms (order >= 2) from k confounders
n_interactions_from_k <- function(k) {

  # if k < 2, there are no interactions
  if (k < 2) return(0L)

  # number of non-linear coefficients: all interactions of order >= 2
  sum(sapply(2:k, function(m) choose(k, m)))
}

# helper: coerce a user-provided B input into a list of length T
# accepts either:
#   - a single 2 x p matrix (recycled to all waves), or
#   - a list of length T, each element a 2 x p matrix
as_B_list <- function(B_in, T) {

  # if user provides a single matrix, recycle it across T waves
  if (is.matrix(B_in)) {

    # basic sanity check
    if (nrow(B_in) != 2)
      stop("B matrix must have exactly 2 rows (X, Y).")

    B_list <- vector("list", T)
    names(B_list) <- paste0("t", 1:T)
    for (t in 1:T) B_list[[t]] <- B_in
    return(B_list)
  }

  # if user provides a list, validate length and contents
  if (is.list(B_in)) {

    if (length(B_in) != T)
      stop("B list must have length T (one 2 x p matrix per wave).")

    for (t in 1:T) {
      if (!is.matrix(B_in[[t]]))
        stop("Each element of B_list must be a 2 x p matrix.")
      if (nrow(B_in[[t]]) != 2)
        stop("Each B_t must have exactly 2 rows (X, Y).")
    }

    # set names if missing
    if (is.null(names(B_in))) names(B_in) <- paste0("t", 1:T)

    return(B_in)
  }

  stop("B_in must be either a 2 x p matrix or a list of length T of 2 x p matrices.")
}

############################################################
##  3. Simulate panel data
############################################################

# this function simulates panel data under a CLPM data-generating process, where 1) 
# 1) we simulate the confounders first for each wave, based on the given Psi covariance matrix
# 2) we compute the variance these confounders induce at each wave given the B matrix for that wave
# 3) we compute how much variance must then be induced by the dynamic process to reach the target covariance 
# 4) however, we want some extra covariance between X and Y at each occasion, meaning that we need to first compute
#    the implied covariance at each wave, and then add this extra covariance to the target covariance at each wave
# 5) calculate the innovations covariance needed to reach this target covariance at each wave
# 6) simulate the panel data panel data for each wave using these innovations
#    then adding the lagged effects from A and the direct confounder effects from B
# 7) repeat for T waves with varying B matrices

simulate_panel_data_int <- function(
    N,                                                         # number of individuals
    T,                                                         # number of waves
    A,                                                         # 2x2 autoregressive/cross-lag matrix
    B_list,                                                    # list of B matrices: B_list[[t]] is 2 x p
    Psi,                                                       # k x k confounder covariance matrix (linear confounders)
    rho_extra,                                                 # extra covariance to add at observed level
    seed = NULL                                                # optional random seed for reproducibility
){

  # helper function to find stationary covariance c given A and S_U
  # given that the innovations are uncorrelated, what is the covariance between X and Y?
  # we need this to be able to add some extra covariance at the observed level on top of the existing covariance
  find_c <- function(A, S_U) {
    
    # given A and S_U, find the stationary covariance c between X and Y
    f <- function(c) {
      # given a candidate covariance c, compute the correlation of the innovations
      # we want this correlation to be 0, so we try to find a c that achieves this

      # target stationary covariance for (X_t, Y_t)
      # where the variance is 1 (as always) and covariance is c
      S_target_c <- matrix(c(1, c,
                             c, 1),
                           nrow = 2, byrow = TRUE)

      # dynamic component variance: S_dyn = S_target - S_U
      S_dyn_c <- S_target_c - S_U

      # ensure symmetry
      S_dyn_c <- (S_dyn_c + t(S_dyn_c)) / 2

      # innovations covariance implied by stationarity:
      # S_dyn = A S_dyn A' + Sigma_e_c
      Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A

      # ensure symmetry
      Sigma_e_c <- (Sigma_e_c + t(Sigma_e_c)) / 2

      # compute correlation of innovations (should match rho = 0)
      v1    <- Sigma_e_c[1, 1]
      v2    <- Sigma_e_c[2, 2]
      cov12 <- Sigma_e_c[1, 2]
      corr_e <- cov12 / sqrt(v1 * v2)

      corr_e - 0   # we want corr(e_x, e_y) = 0
    }

    # root-finding for covariance c between -0.99 and 0.99
    # wrapped in tryCatch so the simulation does not break if no root exists
    out <- tryCatch(
      uniroot(f, interval = c(-0.99, 0.99))$root,
      error = function(e) NA_real_
    )

    out
  }

  # set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # number of confounders (linear)
  k <- ncol(Psi)

  # simulate the confounders from their multivariate normal distribution
  # given a specified covariance matrix Psi, where diag=1 and for now the off-diag=0, giving I, the identity matrix. 
  U <- mvtnorm::rmvnorm(
    n     = N,                                                 # number of individuals
    mean  = rep(0, k),                                         # all confounders have mean 0
    sigma = Psi                                                # covariance matrix of confounders
  )

  # generate all interaction index sets of order >= 2
  int_idx <- if (k >= 2) {
    unlist(
      # for each order m from 2 to k
      # combine k confounders taken m at a time
      lapply(2:k, function(m) combn(1:k, m, simplify = FALSE)),
      # flatten the list
      recursive = FALSE
    )
  } else {
    list()
  }

  # build raw interaction terms
  U_int_raw <- if (length(int_idx) > 0) {
    sapply(int_idx, function(ix) {
      # compute product of confounders in this interaction
      apply(U[, ix, drop = FALSE], 1, prod)
    })
  } else {
    matrix(nrow = N, ncol = 0)
  }

  # orthogonalise interaction terms with respect to their linear parents
  # and normalise them to variance 1
  U_int_orth <- if (length(int_idx) > 0) {
    sapply(seq_along(int_idx), function(j) {

      # extract the parents for this interaction
      parents <- U[, int_idx[[j]], drop = FALSE]

      # raw interaction term
      z_raw <- U_int_raw[, j]

      # residualise the interaction on its parents
      z_orth <- resid(lm(z_raw ~ parents))

      # normalise to variance 1
      sdz <- sd(z_orth)
      if (is.na(sdz) || sdz == 0) z_orth else z_orth / sdz
    })
  } else {
    matrix(nrow = N, ncol = 0)
  }

  # combine linear and non-linear confounders
  U_full <- cbind(U, U_int_orth)

  
  # preparing containers for the variance structure
  S_dyn_list    <- vector("list", T)                           # variance coming from crosslaggs, autoreg and innovations
  Sigma_e_list  <- vector("list", T)                           # innovations covariance
  S_target_list <- vector("list", T)                           # target covariance at observed level
  c_base_vec    <- numeric(T)                                  # the baseline covariance if residual covariance is 0 (meaning all coming from system + confounders)
  c_total_vec   <- numeric(T)                                  # the total covariance at observed level (including extra rho)

  # computing the variance structure at each wave
  # for each wave t from 1 to T
  for (t in 1:T) {

    # get the B matrix for this wave
    B_t <- B_list[[t]]

    # the confounder induced variance covariance at this wave is B_t Psi B_t'
    # note: Psi is implicitly the identity for the full confounder set due to orthogonalisation
    S_U_t <- B_t %*% t(B_t)

    # the base covariance can be found (if the innovations are uncorrelated) by calling find_c
    c_base_t <- find_c(A, S_U_t)

    # store the base covariance in the container
    c_base_vec[t] <- c_base_t

    # now we need to add some extra covariance at the observed level
    c_total_t <- c_base_t + rho_extra

    # and store this in the container
    c_total_vec[t] <- c_total_t

    # S_target can then be specified using our computed total covariance and variance = 1
    S_target_t <- matrix(c(1, c_total_t,
                           c_total_t, 1),
                         nrow = 2, byrow = TRUE)
    
    # and save this in the container
    S_target_list[[t]] <- S_target_t

    # the variance coming from the dynamic process (cross-lag + auto + innovations) is then:
    # target - confounder induced variance
    S_dyn_t <- S_target_t - S_U_t

    # ensure symmetry
    S_dyn_t <- (S_dyn_t + t(S_dyn_t))/2   

    # store in container
    S_dyn_list[[t]] <- S_dyn_t

    # innovations covariance from stationarity:
    # S_dyn = A S_dyn A' + Sigma_e
    Sigma_e_t <- S_dyn_t - t(A) %*% S_dyn_t %*% A

    # ensure symmetry
    Sigma_e_t <- (Sigma_e_t + t(Sigma_e_t))/2

    # store in container
    Sigma_e_list[[t]] <- Sigma_e_t
  }

  # prepare data frame to hold simulated data
  # number of confounders
  p <- ncol(U_full)

  # data frame with N rows and 2*T + p columns
  df <- matrix(NA, nrow = N, ncol = 2*T + p)

  # set column names
  colnames(df) <- c(paste0("x", 1:T),
                    paste0("y", 1:T),
                    paste0("c", 1:p))

  # add confounders to dataframe
  df[, (2*T + 1):(2*T + p)] <- U_full

  # simulate the first wave, which is different because there are no lagged values yet
  Ddyn <- mvtnorm::rmvnorm(
    n     = N,                                                   # number of individuals
    mean  = c(0, 0),                                             # mean 0 for X and Y
    sigma = S_dyn_list[[1]]                                      # variance covariance matrix at wave 1
  ) 

  # add the direct confounder effects
  obs1 <- Ddyn + U_full %*% t(B_list[[1]])

  # store in dataframe
  df[, "x1"] <- obs1[, 1]
  df[, "y1"] <- obs1[, 2]

  # simulate waves 2 to T
  for (t in 2:T) {

    # pull variance covariance matrix for this wave
    Sigma_e_t <- Sigma_e_list[[t]]

    # dynamic process
    Ddyn <- Ddyn %*% t(A) + mvtnorm::rmvnorm(N, sigma = Sigma_e_t)

    # add direct confounder effects
    obs <- Ddyn + U_full %*% t(B_list[[t]])

    # store
    df[, paste0("x", t)] <- obs[, 1]
    df[, paste0("y", t)] <- obs[, 2]
  }

  return(df)
}

############################################################
##  4. Model builders
############################################################

# we want our models to adapt to the number of time points T, and since those models are strings
# we will need to built them using text manipulation

# CLPM model string builder, without confounder adjustment at all
build_clpm <- function(T) {

  # here we build the lines:
  # X_t = X_{t-1} + Y_{t-1}
  # Y_t = X_{t-1} + Y_{t-1}
  regress_block <- paste(

    # for each time point from 2 to T
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )

    # add a line break between each time point
    })), collapse="\n"
  )

  # now we need to add the residual covariances
  # producing X_t ~~ Y_t
  resid_cov <- paste(sprintf("x%d ~~ y%d", 1:T, 1:T), collapse="\n")

  # the residual variances for X_t and Y_t
  resid_vars <- paste(

    # yielding lines like X_t ~~ X_t
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),

    # and Y_t ~~ Y_t
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  # we now need to set the means to 1
  means_block <- paste(

    # produces lines: x1 + x2 + ... + xT ~ 1
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",

    # produces lines: y1 + y2 + ... + yT ~ 1
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  # combine all blocks into one model string
  paste(regress_block, resid_cov, resid_vars, means_block, sep="\n")
}

# same as above, but with direct confounder adjustment added
build_clpm_with_Cs <- function(T, k) {

  # creates the line c1 + c2 + ... + ck
  C_names <- paste0("c", 1:k, collapse=" + ")

  # autoregressive and cross-lagged paths, but also confounders added
  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # produces: X_t ~ X_{t-1} + Y_{t-1} + c1 + c2 + ... + ck
        sprintf("x%d ~ x%d + y%d + %s", t, t-1, t-1, C_names),

        # produces: Y_t ~ X_{t-1} + Y_{t-1} + c1 + c2 + ... + ck
        sprintf("y%d ~ x%d + y%d + %s", t, t-1, t-1, C_names)
      )
    })), collapse="\n"
  )

  # from here the function behaves the same as above
  resid_cov <- paste(sprintf("x%d ~~ y%d", 1:T, 1:T), collapse="\n")

  resid_vars <- paste(
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  means_block <- paste(
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  paste(regress_block, resid_cov, resid_vars, means_block, sep="\n")
}

# same as above, but with indirect confounder adjustment via random intercepts
build_riclpm <- function(T) {

  # here we create the random intercepts
  ri_block <- paste0(

    # produces lines like rix =~ 1*x1 + 1*x2 + ... + 1*xT
    "rix =~ ", paste(sprintf("1*x%d", 1:T), collapse=" + "), "\n",

    # produces lines like riy =~ 1*y1 + 1*y2 + ... + 1*yT
    "riy =~ ", paste(sprintf("1*y%d", 1:T), collapse=" + "), "\n",

    # since this is allways the same, we directly add the variances and covariance of the random intercepts
    "rix ~~ rix\n riy ~~ riy\n rix ~~ riy\n"
  )

  # here we fix the residual variances to zero
  resid_fix <- paste0(

    # produces lines like x1 ~~ 0*x1 + 0*x2 + ... + 0*xT
    paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse="; "), "\n",

    # and y1 ~~ 0*y1 + 0*y2 + ... + 0*yT
    paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the within-person latent variables for X_t and Y_t
  within_lat <- paste0(

    # produces lines like wx1 =~ 1*x1, wx2 =~ 1*x2, ..., wxT =~ 1*xT
    paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse="; "), "\n",

    # and wy1 =~ 1*y1, wy2 =~ 1*y2, ..., wyT =~ 1*yT
    paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the orthogonality constraints: i.e. stable traits are uncorrelated with within-person fluctuations
  orth <- paste0(
    "rix ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "rix ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n"
  )

  # here we create the within-person variances
  within_var <- paste0(

    # creates lines like wx1 ~~ wx1, wx2 ~~ wx2, ..., wxT ~~ wxT
    paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n",

    # and wy1 ~~ wy1, wy2 ~~ wy2, ..., wyT ~~ wyT
    paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the within-person covariances
  within_cov <- paste0(

    # creates lines like wx1 ~~ wy1, wx2 ~~ wy2, ..., wxT ~~ wyT
    paste(sprintf("wy%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the autoregressive and cross-lagged paths
  regress <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}: wx_t ~ wx_{t-1} + wy_{t-1}
        sprintf("wx%d ~ wx%d + wy%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}: wy_t ~ wx_{t-1} + wy_{t-1}
        sprintf("wy%d ~ wx%d + wy%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  # here we create the means
  means <- paste0(

    # produces lines like x1 ~ mx*1, y1 ~ my*1
    paste(paste0("x",1:T), collapse=" + "), " ~ mx*1\n",

    # produces lines like x1 ~ mx*1, y1 ~ my*1
    paste(paste0("y",1:T), collapse=" + "), " ~ my*1\n"
  )

  # finally, we put it all together
  paste(ri_block, resid_fix, within_lat, orth,
        within_var, within_cov, regress, means, sep="\n")
}

# now we build the DPM model string builder
build_dpm <- function(T) {

  # define the accumulating factors FX 
  FX_block <- paste0(

    # produces line FX =~ 1*x1 + 1*x2 + ... + 1*xT
    "FX =~ ", paste(sprintf("1*x%d", 2:T), collapse=" + "), "\n"
  )

  # define the accumulating factors FY
  FY_block <- paste0(

    # produces line FY =~ 1*y1 + 1*y2 + ... + 1*yT
    "FY =~ ", paste(sprintf("1*y%d", 2:T), collapse=" + "), "\n"
  )

  # define the residual covariances between FX and x1, and FY and y1
  fx_cov_block <- "FX ~~ x1 + y1\n"
  fy_cov_block <- "FY ~~ x1 + y1\n"

  # define the autoregressive and cross-lagged paths
  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  # define the residual covariances between X_t and Y_t
  resid_cov_block <- paste(

    # produces lines like X_t ~~ Y_t
    sprintf("x%d ~~ y%d", 1:T, 1:T),
    collapse="\n"
  )

  # define the latent covariances between FX and FY
  latent_cov_block <- paste(
    "FX ~~ FX",
    "FY ~~ FY",
    "FX ~~ FY",
    sep="\n"
  )

  # define the residual variances
  resid_var_block <- paste(

    # produces lines like X_t ~~ X_t
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),

    # produces lines like Y_t ~~ Y_t
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  # define the means
  means_block <- paste(

    # produces lines like x1 ~ 1, y1 ~ 1
    paste(sprintf("x%d", 1:T), collapse=" + "), "~ 1\n",

    # produces lines like x1 ~ 1, y1 ~ 1
    paste(sprintf("y%d", 1:T), collapse=" + "), "~ 1\n"
  )

  # finally, we put it all together
  paste(
    FX_block,
    FY_block,
    fx_cov_block,
    fy_cov_block,
    regress_block,
    resid_cov_block,
    latent_cov_block,
    resid_var_block,
    means_block,
    sep="\n"
  )
}

############################################################
##  5. Residualizer 
############################################################

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

############################################################
##  6. SAFE FITTING HELPERS (capture error messages)
############################################################

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

# same as above but for CLPM with confounders
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

# same as above but for CLPM with XGB residualisation
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

############################################################
##  6.2 XGB: One-time CV tuner (run once per study)
############################################################

# in a true application of XBG residualisation, one would want to tune the hyperparameters
# however, this is computationally very costly. 
# therefore, we provide a function that can be run once per study design to find the optimal hyperparameters
# and recycle these for all simulations in that design, recreating how tuning would happen in practice

tune_xgb_once <- function(df,
                          k,
                          x_prefix = "x",
                          y_prefix = "y",
                          c_prefix = "c",
                          nfold = 5,
                          nrounds_max = 6000,                 # (CHANGED) allow more rounds for smaller eta
                          early_stopping_rounds = 100,        # (CHANGED) reduce "too-early stop"
                          max_grid = 600,                     # (NEW) slightly larger search budget
                          seed = 1) {                         # (NEW) reproducible grid + folds

  # convert to data frame
  df <- as.data.frame(df)

  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)

  # only linear confounders are observed
  c_cols <- paste0(c_prefix, 1:k)

  # confounder matrix (standardised for XGB stability)
  C <- scale(as.matrix(df[c_cols]))

  # (NEW) helper: tune on a stacked target with grouped folds by person
  tune_one <- function(y_stack, C_stack, n_person) {

    dtrain <- xgboost::xgb.DMatrix(data = C_stack, label = y_stack)

    # (NEW) grouped folds (avoid leakage when the same person's C row is repeated across waves)
    person_id <- rep(seq_len(n_person), times = nrow(C_stack) / n_person)

    set.seed(seed)
    fold_id_person <- sample(rep(1:nfold, length.out = n_person))

    folds <- lapply(1:nfold, function(f) which(fold_id_person[person_id] == f))

    # (CHANGED) expanded grid (still “small but useful”, just a bit richer)
    grid <- expand.grid(
      max_depth        = c(2:10),
      min_child_weight = c(1, 2, 5, 10, 20),
      eta              = c(0.003, 0.005, 0.01, 0.02, 0.05, 0.10),
      subsample        = c(0.5, 0.6, 0.8, 1.0),
      colsample_bytree = c(0.5, 0.6, 0.8, 1.0),
      gamma            = c(0, 0.1, 0.5, 1, 2, 5),
      lambda           = c(0.5, 1, 2, 5),
      alpha            = c(0, 0.001, 0.01, 0.1),
      KEEP.OUT.ATTRS   = FALSE
    )

    # trim the grid to keep tuning feasible
    if (nrow(grid) > max_grid) {
      set.seed(seed)
      grid <- grid[sample.int(nrow(grid), max_grid), , drop = FALSE]
    }

    best_rmse    <- Inf
    best_params  <- NULL
    best_nrounds <- NULL

    # create a progress bar for the CV grid
    # (this one shows elapsed time + estimated time left)
    pb <- pbapply::timerProgressBar(min = 0, max = nrow(grid), style = 3)

    for (i in seq_len(nrow(grid))) {

      # update progress bar
      pbapply::setTimerProgressBar(pb, i)

      params <- as.list(grid[i, ])

      cv <- xgboost::xgb.cv(
        data = dtrain,
        nrounds = nrounds_max,
        folds = folds,                          # (NEW) grouped folds to avoid leakage
        early_stopping_rounds = early_stopping_rounds,
        verbose = 0,
        params = c(
          list(
            objective = "reg:squarederror",
            eval_metric = "rmse",
            booster = "gbtree",
            tree_method = "hist",

            # let xgboost use multiple threads during tuning
            nthread = max(1, parallel::detectCores() - 1)
          ),
          params
        )
      )

      rmse <- cv$evaluation_log$test_rmse_mean[cv$best_iteration]

      if (rmse < best_rmse) {
        best_rmse    <- rmse
        best_params  <- params
        best_nrounds <- cv$best_iteration
      }
    }

    # close the progress bar
    close(pb)

    list(
      params  = best_params,
      nrounds = best_nrounds,
      rmse    = best_rmse
    )
  }

  # message so it is clear why the simulation progress bar has not started yet
  cat("\nTuning the XGB model (one-time CV)...\n")

  # (NEW) tune X and Y separately
  # stack all X waves to tune on more rows
  y_stack_x <- unlist(df[x_cols])
  C_stack_x <- C[rep(seq_len(nrow(C)), times = length(x_cols)), , drop = FALSE]

  cat("\n - Tuning X residualiser...\n")
  tuned_x <- tune_one(y_stack_x, C_stack_x, n_person = nrow(C))

  # stack all Y waves to tune on more rows
  y_stack_y <- unlist(df[y_cols])
  C_stack_y <- C[rep(seq_len(nrow(C)), times = length(y_cols)), , drop = FALSE]

  cat("\n - Tuning Y residualiser...\n")
  tuned_y <- tune_one(y_stack_y, C_stack_y, n_person = nrow(C))

  # newline so the next progress bar starts on a clean line
  cat("\n")

  list(
    params_X  = tuned_x$params,
    nrounds_X = tuned_x$nrounds,
    rmse_X    = tuned_x$rmse,

    params_Y  = tuned_y$params,
    nrounds_Y = tuned_y$nrounds,
    rmse_Y    = tuned_y$rmse
  )
}

############################################################
##  7. Extract lagged parameters + rho
############################################################

# since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive, cross-lagged, and residual correlation (rho) are located

# lagged parameters extractor
# (CHANGED) now returns: est, p-value, CI lower, CI upper for each lagged parameter
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

############################################################
##  8. Wrapper for one replication
############################################################

# so now we need to obtain a function that runs all other functions for one replication of the study
# this function should contain all the arguments that are defined in the previous sections

# (CHANGED) Instead of sampling/generating B, the user passes B_scenarios:
#   - named list: B_scenarios[["constant"]] = (matrix or list length T)
#   - or unnamed list aligned with `scenarios` order
run_one_rep_study2 <- function(
    rep_id,                                                 # replication index (set by outer loop)
    N,                                                      # sample size
    T,                                                      # number of waves
    k,                                                      # number of linear confounder parts
    scenarios,                                              # character vector: c("constant","stepwise",...)
    B_scenarios,                                            # (NEW) user-supplied B matrices per scenario
    A,                                                      # 2×2 autoregressive + cross-lag matrix
    Psi,                                                    # k×k confounder covariance
    rho_extra,                                              # extra covariance added to X,Y each wave
    models_to_run,                                          # e.g. c("clpm","riclpm","dpm","lbca","adj","xgb")
    base_seed = 1234,                                       # base seed
    ci_level = 0.95                                         # CI level for extracted parameters
){

  # set seed for this replication, use rep_id so it varies
  set.seed(base_seed + rep_id)

  # expected number of confounders (linear + interactions)
  k_non <- n_interactions_from_k(k)
  p_exp <- k + k_non

  # prepare output list for each scenario
  out_list <- vector("list", length(scenarios))

  # loop over the scenarios
  for (j in seq_along(scenarios)) {

    scen <- scenarios[j]

    # select scenario-specific B input
    # supports:
    #   - named list (preferred): B_scenarios[[scen]]
    #   - unnamed list aligned with scenarios: B_scenarios[[j]]
    B_in <- NULL
    if (is.list(B_scenarios) && !is.null(names(B_scenarios)) && scen %in% names(B_scenarios)) {
      B_in <- B_scenarios[[scen]]
    } else if (is.list(B_scenarios) && is.null(names(B_scenarios)) && length(B_scenarios) == length(scenarios)) {
      B_in <- B_scenarios[[j]]
    } else if (is.list(B_scenarios) && length(B_scenarios) == 1L) {
      # if only one B provided, recycle across scenarios
      B_in <- B_scenarios[[1]]
    } else {
      stop("Could not map B_scenarios to scenarios. Provide a named list with scenario names, or an unnamed list aligned with `scenarios`.")
    }

    # coerce to list-of-length-T
    B_list <- as_B_list(B_in, T)

    # validate dimensions vs k and Psi
    if (!all(dim(Psi) == c(k, k)))
      stop("Mismatch: Psi must be k x k with k = ", k, ".")
    for (t in 1:T) {
      if (ncol(B_list[[t]]) != p_exp)
        stop("Mismatch: B_list[[", t, "]] has p = ", ncol(B_list[[t]]),
             " columns, but expected p = k + #interactions = ", p_exp, ".")
    }

    # extract mean betas
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec

    # try to simulate panel data
    df <- tryCatch(
      simulate_panel_data_int(
        N         = N,
        T         = T,
        A         = A,
        B_list    = B_list,
        Psi       = Psi,
        rho_extra = rho_extra
      ),
      error = function(e) NULL
    )

    if (is.null(df)) {

      # simulation failed return NA
      out_list[[j]] <- data.frame(
        run      = rep(rep_id, T),
        occasion = 1:T,
        scenario = scen,
        beta     = beta_vec,
        beta_X   = beta_X_vec,
        beta_Y   = beta_Y_vec,

        # cross-lag XY (est + p + CI)
        estXY_CLPM      = NA, pXY_CLPM      = NA, ciL_XY_CLPM      = NA, ciU_XY_CLPM      = NA,
        estXY_RI_CLPM   = NA, pXY_RI_CLPM   = NA, ciL_XY_RI_CLPM   = NA, ciU_XY_RI_CLPM   = NA,
        estXY_DPM       = NA, pXY_DPM       = NA, ciL_XY_DPM       = NA, ciU_XY_DPM       = NA,
        estXY_CLPM_Adj  = NA, pXY_CLPM_Adj  = NA, ciL_XY_CLPM_Adj  = NA, ciU_XY_CLPM_Adj  = NA,
        estXY_CLPM_LBCA = NA, pXY_CLPM_LBCA = NA, ciL_XY_CLPM_LBCA = NA, ciU_XY_CLPM_LBCA = NA,
        estXY_CLPM_XGB  = NA, pXY_CLPM_XGB  = NA, ciL_XY_CLPM_XGB  = NA, ciU_XY_CLPM_XGB  = NA,

        # cross-lag YX (est + p + CI)
        estYX_CLPM      = NA, pYX_CLPM      = NA, ciL_YX_CLPM      = NA, ciU_YX_CLPM      = NA,
        estYX_RI_CLPM   = NA, pYX_RI_CLPM   = NA, ciL_YX_RI_CLPM   = NA, ciU_YX_RI_CLPM   = NA,
        estYX_DPM       = NA, pYX_DPM       = NA, ciL_YX_DPM       = NA, ciU_YX_DPM       = NA,
        estYX_CLPM_Adj  = NA, pYX_CLPM_Adj  = NA, ciL_YX_CLPM_Adj  = NA, ciU_YX_CLPM_Adj  = NA,
        estYX_CLPM_LBCA = NA, pYX_CLPM_LBCA = NA, ciL_YX_CLPM_LBCA = NA, ciU_YX_CLPM_LBCA = NA,
        estYX_CLPM_XGB  = NA, pYX_CLPM_XGB  = NA, ciL_YX_CLPM_XGB  = NA, ciU_YX_CLPM_XGB  = NA,

        # AR X (est + p + CI)
        estA_CLPM      = NA, pA_CLPM      = NA, ciL_A_CLPM      = NA, ciU_A_CLPM      = NA,
        estA_RI_CLPM   = NA, pA_RI_CLPM   = NA, ciL_A_RI_CLPM   = NA, ciU_A_RI_CLPM   = NA,
        estA_DPM       = NA, pA_DPM       = NA, ciL_A_DPM       = NA, ciU_A_DPM       = NA,
        estA_CLPM_Adj  = NA, pA_CLPM_Adj  = NA, ciL_A_CLPM_Adj  = NA, ciU_A_CLPM_Adj  = NA,
        estA_CLPM_LBCA = NA, pA_CLPM_LBCA = NA, ciL_A_CLPM_LBCA = NA, ciU_A_CLPM_LBCA = NA,
        estA_CLPM_XGB  = NA, pA_CLPM_XGB  = NA, ciL_A_CLPM_XGB  = NA, ciU_A_CLPM_XGB  = NA,

        # AR Y (est + p + CI)
        estAY_CLPM      = NA, pAY_CLPM      = NA, ciL_AY_CLPM      = NA, ciU_AY_CLPM      = NA,
        estAY_RI_CLPM   = NA, pAY_RI_CLPM   = NA, ciL_AY_RI_CLPM   = NA, ciU_AY_RI_CLPM   = NA,
        estAY_DPM       = NA, pAY_DPM       = NA, ciL_AY_DPM       = NA, ciU_AY_DPM       = NA,
        estAY_CLPM_Adj  = NA, pAY_CLPM_Adj  = NA, ciL_AY_CLPM_Adj  = NA, ciU_AY_CLPM_Adj  = NA,
        estAY_CLPM_LBCA = NA, pAY_CLPM_LBCA = NA, ciL_AY_CLPM_LBCA = NA, ciU_AY_CLPM_LBCA = NA,
        estAY_CLPM_XGB  = NA, pAY_CLPM_XGB  = NA, ciL_AY_CLPM_XGB  = NA, ciU_AY_CLPM_XGB  = NA,

        # residual correlation
        # (NEW) now also p-value and CI for all methods
        estRho_CLPM      = NA, pRho_CLPM      = NA, ciL_Rho_CLPM      = NA, ciU_Rho_CLPM      = NA,
        estRho_RI_CLPM   = NA, pRho_RI_CLPM   = NA, ciL_Rho_RI_CLPM   = NA, ciU_Rho_RI_CLPM   = NA,
        estRho_DPM       = NA, pRho_DPM       = NA, ciL_Rho_DPM       = NA, ciU_Rho_DPM       = NA,
        estRho_CLPM_Adj  = NA, pRho_CLPM_Adj  = NA, ciL_Rho_CLPM_Adj  = NA, ciU_Rho_CLPM_Adj  = NA,
        estRho_CLPM_LBCA = NA, pRho_CLPM_LBCA = NA, ciL_Rho_CLPM_LBCA = NA, ciU_Rho_CLPM_LBCA = NA,
        estRho_CLPM_XGB  = NA, pRho_CLPM_XGB  = NA, ciL_Rho_CLPM_XGB  = NA, ciU_Rho_CLPM_XGB  = NA,

        fail_CLPM      = TRUE,
        fail_RI_CLPM   = TRUE,
        fail_DPM       = TRUE,
        fail_CLPM_Adj  = TRUE,
        fail_CLPM_LBCA = TRUE,
        fail_CLPM_XGB  = TRUE,

        err_CLPM      = "sim failed",
        err_RI_CLPM   = "sim failed",
        err_DPM       = "sim failed",
        err_CLPM_Adj  = "sim failed",
        err_CLPM_LBCA = "sim failed",
        err_CLPM_XGB  = "sim failed",

        is_na_run = 1L
      )

      next
    }

    # build model strings
    model_clpm         <- build_clpm(T)
    model_riclpm       <- build_riclpm(T)
    model_dpm          <- build_dpm(T)
    model_clpm_with_Cs <- build_clpm_with_Cs(T, k)

    # fit the models safely
    res_clpm <- if ("clpm"   %in% models_to_run) safe_fit_clpm(model_clpm, df) else list(fit=NULL, err=NA)
    res_ric  <- if ("riclpm" %in% models_to_run) safe_fit_riclpm(model_riclpm, df) else list(fit=NULL, err=NA)
    res_dpm0 <- if ("dpm"    %in% models_to_run) safe_fit_dpm(model_dpm, df) else list(fit=NULL, err=NA)
    res_adj  <- if ("adj"    %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df) else list(fit=NULL, err=NA)

    # LBCA + XGB should only use the linear confounders (c1..ck)
    res_lbca <- if ("lbca"   %in% models_to_run) safe_fit_clpm_resid(model_clpm, df, k) else list(fit=NULL, err=NA)
    res_xgb  <- if ("xgb"    %in% models_to_run) safe_fit_clpm_xgb(model_clpm, df, k) else list(fit=NULL, err=NA)

    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit
    fit_xgb      <- res_xgb$fit

    # extract lagged parameters (NOW includes p-values + CI)
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",    ci_level = ci_level)
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm",  ci_level = ci_level)
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm",     ci_level = ci_level)
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm",    ci_level = ci_level)
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm",    ci_level = ci_level)
    lag_xgb  <- extract_lagged_parameters(fit_xgb,       T, "clpm",    ci_level = ci_level)

    # residual correlations
    # (NEW) now returns a data.frame with est, p, CI for all methods
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm",   ci_level = ci_level)
    rho_ric  <- extract_rho_vec(fit_ric,      T, "riclpm", ci_level = ci_level)
    rho_dpm  <- extract_rho_vec(fit_dpm0,     T, "dpm",    ci_level = ci_level)
    rho_adj  <- extract_rho_vec(fit_adj,      T, "clpm",   ci_level = ci_level)
    rho_lbca <- extract_rho_vec(fit_lbca,     T, "clpm",   ci_level = ci_level)
    rho_xgb  <- extract_rho_vec(fit_xgb,      T, "clpm",   ci_level = ci_level)

    # assemble output row
    out_list[[j]] <- data.frame(

      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      beta     = beta_vec,
      beta_X   = beta_X_vec,
      beta_Y   = beta_Y_vec,

      # cross-lag XY (X → Y)
      estXY_CLPM      = c(NA, lag_raw$xy$est),
      pXY_CLPM        = c(NA, lag_raw$xy$p),
      ciL_XY_CLPM     = c(NA, lag_raw$xy$ci.lower),
      ciU_XY_CLPM     = c(NA, lag_raw$xy$ci.upper),

      estXY_RI_CLPM   = c(NA, lag_ric$xy$est),
      pXY_RI_CLPM     = c(NA, lag_ric$xy$p),
      ciL_XY_RI_CLPM  = c(NA, lag_ric$xy$ci.lower),
      ciU_XY_RI_CLPM  = c(NA, lag_ric$xy$ci.upper),

      estXY_DPM       = c(NA, lag_dpm0$xy$est),
      pXY_DPM         = c(NA, lag_dpm0$xy$p),
      ciL_XY_DPM      = c(NA, lag_dpm0$xy$ci.lower),
      ciU_XY_DPM      = c(NA, lag_dpm0$xy$ci.upper),

      estXY_CLPM_Adj  = c(NA, lag_adj$xy$est),
      pXY_CLPM_Adj    = c(NA, lag_adj$xy$p),
      ciL_XY_CLPM_Adj = c(NA, lag_adj$xy$ci.lower),
      ciU_XY_CLPM_Adj = c(NA, lag_adj$xy$ci.upper),

      estXY_CLPM_LBCA  = c(NA, lag_lbca$xy$est),
      pXY_CLPM_LBCA    = c(NA, lag_lbca$xy$p),
      ciL_XY_CLPM_LBCA = c(NA, lag_lbca$xy$ci.lower),
      ciU_XY_CLPM_LBCA = c(NA, lag_lbca$xy$ci.upper),

      estXY_CLPM_XGB  = c(NA, lag_xgb$xy$est),
      pXY_CLPM_XGB    = c(NA, lag_xgb$xy$p),
      ciL_XY_CLPM_XGB = c(NA, lag_xgb$xy$ci.lower),
      ciU_XY_CLPM_XGB = c(NA, lag_xgb$xy$ci.upper),

      # cross-lag YX (Y → X)
      estYX_CLPM      = c(NA, lag_raw$yx$est),
      pYX_CLPM        = c(NA, lag_raw$yx$p),
      ciL_YX_CLPM     = c(NA, lag_raw$yx$ci.lower),
      ciU_YX_CLPM     = c(NA, lag_raw$yx$ci.upper),

      estYX_RI_CLPM   = c(NA, lag_ric$yx$est),
      pYX_RI_CLPM     = c(NA, lag_ric$yx$p),
      ciL_YX_RI_CLPM  = c(NA, lag_ric$yx$ci.lower),
      ciU_YX_RI_CLPM  = c(NA, lag_ric$yx$ci.upper),

      estYX_DPM       = c(NA, lag_dpm0$yx$est),
      pYX_DPM         = c(NA, lag_dpm0$yx$p),
      ciL_YX_DPM      = c(NA, lag_dpm0$yx$ci.lower),
      ciU_YX_DPM      = c(NA, lag_dpm0$yx$ci.upper),

      estYX_CLPM_Adj  = c(NA, lag_adj$yx$est),
      pYX_CLPM_Adj    = c(NA, lag_adj$yx$p),
      ciL_YX_CLPM_Adj = c(NA, lag_adj$yx$ci.lower),
      ciU_YX_CLPM_Adj = c(NA, lag_adj$yx$ci.upper),

      estYX_CLPM_LBCA  = c(NA, lag_lbca$yx$est),
      pYX_CLPM_LBCA    = c(NA, lag_lbca$yx$p),
      ciL_YX_CLPM_LBCA = c(NA, lag_lbca$yx$ci.lower),
      ciU_YX_CLPM_LBCA = c(NA, lag_lbca$yx$ci.upper),

      estYX_CLPM_XGB  = c(NA, lag_xgb$yx$est),
      pYX_CLPM_XGB    = c(NA, lag_xgb$yx$p),
      ciL_YX_CLPM_XGB = c(NA, lag_xgb$yx$ci.lower),
      ciU_YX_CLPM_XGB = c(NA, lag_xgb$yx$ci.upper),

      # autoregressive X
      estA_CLPM      = c(NA, lag_raw$ar_x$est),
      pA_CLPM        = c(NA, lag_raw$ar_x$p),
      ciL_A_CLPM     = c(NA, lag_raw$ar_x$ci.lower),
      ciU_A_CLPM     = c(NA, lag_raw$ar_x$ci.upper),

      estA_RI_CLPM   = c(NA, lag_ric$ar_x$est),
      pA_RI_CLPM     = c(NA, lag_ric$ar_x$p),
      ciL_A_RI_CLPM  = c(NA, lag_ric$ar_x$ci.lower),
      ciU_A_RI_CLPM  = c(NA, lag_ric$ar_x$ci.upper),

      estA_DPM       = c(NA, lag_dpm0$ar_x$est),
      pA_DPM         = c(NA, lag_dpm0$ar_x$p),
      ciL_A_DPM      = c(NA, lag_dpm0$ar_x$ci.lower),
      ciU_A_DPM      = c(NA, lag_dpm0$ar_x$ci.upper),

      estA_CLPM_Adj  = c(NA, lag_adj$ar_x$est),
      pA_CLPM_Adj    = c(NA, lag_adj$ar_x$p),
      ciL_A_CLPM_Adj = c(NA, lag_adj$ar_x$ci.lower),
      ciU_A_CLPM_Adj = c(NA, lag_adj$ar_x$ci.upper),

      estA_CLPM_LBCA  = c(NA, lag_lbca$ar_x$est),
      pA_CLPM_LBCA    = c(NA, lag_lbca$ar_x$p),
      ciL_A_CLPM_LBCA = c(NA, lag_lbca$ar_x$ci.lower),
      ciU_A_CLPM_LBCA = c(NA, lag_lbca$ar_x$ci.upper),

      estA_CLPM_XGB  = c(NA, lag_xgb$ar_x$est),
      pA_CLPM_XGB    = c(NA, lag_xgb$ar_x$p),
      ciL_A_CLPM_XGB = c(NA, lag_xgb$ar_x$ci.lower),
      ciU_A_CLPM_XGB = c(NA, lag_xgb$ar_x$ci.upper),

      # autoregressive Y
      estAY_CLPM      = c(NA, lag_raw$ar_y$est),
      pAY_CLPM        = c(NA, lag_raw$ar_y$p),
      ciL_AY_CLPM     = c(NA, lag_raw$ar_y$ci.lower),
      ciU_AY_CLPM     = c(NA, lag_raw$ar_y$ci.upper),

      estAY_RI_CLPM   = c(NA, lag_ric$ar_y$est),
      pAY_RI_CLPM     = c(NA, lag_ric$ar_y$p),
      ciL_AY_RI_CLPM  = c(NA, lag_ric$ar_y$ci.lower),
      ciU_AY_RI_CLPM  = c(NA, lag_ric$ar_y$ci.upper),

      estAY_DPM       = c(NA, lag_dpm0$ar_y$est),
      pAY_DPM         = c(NA, lag_dpm0$ar_y$p),
      ciL_AY_DPM      = c(NA, lag_dpm0$ar_y$ci.lower),
      ciU_AY_DPM      = c(NA, lag_dpm0$ar_y$ci.upper),

      estAY_CLPM_Adj  = c(NA, lag_adj$ar_y$est),
      pAY_CLPM_Adj    = c(NA, lag_adj$ar_y$p),
      ciL_AY_CLPM_Adj = c(NA, lag_adj$ar_y$ci.lower),
      ciU_AY_CLPM_Adj = c(NA, lag_adj$ar_y$ci.upper),

      estAY_CLPM_LBCA  = c(NA, lag_lbca$ar_y$est),
      pAY_CLPM_LBCA    = c(NA, lag_lbca$ar_y$p),
      ciL_AY_CLPM_LBCA = c(NA, lag_lbca$ar_y$ci.lower),
      ciU_AY_CLPM_LBCA = c(NA, lag_lbca$ar_y$ci.upper),

      estAY_CLPM_XGB  = c(NA, lag_xgb$ar_y$est),
      pAY_CLPM_XGB    = c(NA, lag_xgb$ar_y$p),
      ciL_AY_CLPM_XGB = c(NA, lag_xgb$ar_y$ci.lower),
      ciU_AY_CLPM_XGB = c(NA, lag_xgb$ar_y$ci.upper),

      # residual correlation
      # (NEW) now output est + p + CI for all methods
      estRho_CLPM      = rho_clpm$est,
      pRho_CLPM        = rho_clpm$p,
      ciL_Rho_CLPM     = rho_clpm$ci.lower,
      ciU_Rho_CLPM     = rho_clpm$ci.upper,

      estRho_RI_CLPM   = rho_ric$est,
      pRho_RI_CLPM     = rho_ric$p,
      ciL_Rho_RI_CLPM  = rho_ric$ci.lower,
      ciU_Rho_RI_CLPM  = rho_ric$ci.upper,

      estRho_DPM       = rho_dpm$est,
      pRho_DPM         = rho_dpm$p,
      ciL_Rho_DPM      = rho_dpm$ci.lower,
      ciU_Rho_DPM      = rho_dpm$ci.upper,

      estRho_CLPM_Adj  = rho_adj$est,
      pRho_CLPM_Adj    = rho_adj$p,
      ciL_Rho_CLPM_Adj = rho_adj$ci.lower,
      ciU_Rho_CLPM_Adj = rho_adj$ci.upper,

      estRho_CLPM_LBCA  = rho_lbca$est,
      pRho_CLPM_LBCA    = rho_lbca$p,
      ciL_Rho_CLPM_LBCA = rho_lbca$ci.lower,
      ciU_Rho_CLPM_LBCA = rho_lbca$ci.upper,

      estRho_CLPM_XGB  = rho_xgb$est,
      pRho_CLPM_XGB    = rho_xgb$p,
      ciL_Rho_CLPM_XGB = rho_xgb$ci.lower,
      ciU_Rho_CLPM_XGB = rho_xgb$ci.upper,

      # failure indicators
      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),
      fail_CLPM_XGB  = is.null(fit_xgb),

      # error messages
      err_CLPM      = rep(res_clpm$err,   T),
      err_RI_CLPM   = rep(res_ric$err,    T),
      err_DPM       = rep(res_dpm0$err,   T),
      err_CLPM_Adj  = rep(res_adj$err,    T),
      err_CLPM_LBCA = rep(res_lbca$err,   T),
      err_CLPM_XGB  = rep(res_xgb$err,    T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy$est, lag_ric$xy$est, lag_dpm0$xy$est,
        lag_adj$xy$est, lag_lbca$xy$est, lag_xgb$xy$est
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}

############################################################
## 9. Main simulation function — PARALLEL
############################################################

run_simulation_study2 <- function(
    reps,                                                                    # number of replications
    N,                                                                       # sample size
    T,                                                                       # number of waves
    k,                                                                       # number of linear confounders
    scenarios,                                                               # e.g., c("constant","stepwise")
    B_scenarios,                                                             # (NEW) user-supplied B matrices per scenario
    A,                                                                       # 2×2 AR + cross-lag matrix
    Psi,                                                                     # k×k confounder covariance
    rho_extra,                                                               # extra covariance added to observations
    models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca","xgb")
    cores = NULL,                                                            # default is detectCores()/2
    base_seed = 1234,                                                        # master seed for reproducible reps
    ci_level = 0.95                                                          # CI level for extracted parameters
) {

  # if the number of cores is not specified, detect and use half of available cores
  if (is.null(cores)) {

    # detect and use half of available cores
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # XGB: One-time CV tuning (run once per study)
  if ("xgb" %in% models_to_run) {

    # choose a pilot scenario for tuning
    scen_tune <- scenarios[1]

    # set seed for pilot so it is reproducible
    set.seed(base_seed + 1)

    # pull pilot B from user-supplied scenarios
    B_in_pilot <- NULL
    if (is.list(B_scenarios) && !is.null(names(B_scenarios)) && scen_tune %in% names(B_scenarios)) {
      B_in_pilot <- B_scenarios[[scen_tune]]
    } else if (is.list(B_scenarios) && is.null(names(B_scenarios)) && length(B_scenarios) >= 1L) {
      B_in_pilot <- B_scenarios[[1]]
    } else {
      stop("For xgb tuning, could not map B_scenarios to the first scenario. Provide a named list or an aligned list.")
    }

    # coerce to list-of-length-T
    B_list_pilot <- as_B_list(B_in_pilot, T)

    # simulate pilot data
    df_pilot <- simulate_panel_data_int(
      N         = N,
      T         = T,
      A         = A,
      B_list    = B_list_pilot,
      Psi       = Psi,
      rho_extra = rho_extra
    )

    # tune and store globally so residualise_panel_xgb can access it
    # (CHANGED) now returns separate tuning for X and Y
    XGB_TUNED <- tune_xgb_once(df_pilot, k)
    assign("XGB_TUNED", XGB_TUNED, envir = .GlobalEnv)

  } else {

    # if we do not run xgb, set the tuning object to NULL
    XGB_TUNED <- NULL
    assign("XGB_TUNED", XGB_TUNED, envir = .GlobalEnv)
  }

  # if cores is 1, run sequentially without parallelization
  if (cores == 1L) {

    # run sequentially
    results_list <- lapply(
      X = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study2(
          rep_id        = rep_id,
          N             = N,
          T             = T,
          k             = k,
          scenarios     = scenarios,
          B_scenarios   = B_scenarios,
          A             = A,
          Psi           = Psi,
          rho_extra     = rho_extra,
          models_to_run = models_to_run,
          base_seed     = base_seed,
          ci_level      = ci_level
        )
      }
    )
    return(dplyr::bind_rows(results_list))
  }

  # make the cluster
  cl <- parallel::makeCluster(cores)

  # load required packages on each worker
  parallel::clusterEvalQ(cl, {
    library(lavaan)
    library(mvtnorm)
    library(xgboost)
    NULL
  })

  # export all necessary functions and variables to the cluster
  parallel::clusterExport(
    cl,
    c(
      # B helpers + simulation
      "n_interactions_from_k",
      "as_B_list",
      "simulate_panel_data_int",

      # model builders
      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",

      # residualisers
      "residualise_panel_linearC",
      "residualise_panel_xgb",     # (CHANGED) now expects separate X/Y tuning in XGB_TUNED

      # safe fitters
      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "safe_fit_clpm_xgb",

      # extractors
      "extract_lagged_parameters",
      "extract_rho_vec",

      # wrapper
      "run_one_rep_study2",

      # arguments
      "N","T","k","scenarios","B_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed","ci_level",

      # tuned xgb settings
      "XGB_TUNED"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X  = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study2(
        rep_id        = rep_id,
        N             = N,
        T             = T,
        k             = k,
        scenarios     = scenarios,
        B_scenarios   = B_scenarios,
        A             = A,
        Psi           = Psi,
        rho_extra     = rho_extra,
        models_to_run = models_to_run,
        base_seed     = base_seed,
        ci_level      = ci_level
      )
    }
  )

  # stop the cluster
  parallel::stopCluster(cl)

  # return the results
  dplyr::bind_rows(results_list)
}
