############################################################
# Simulation Study 1 Functions
############################################################
# 
# The goal of this script is to provide evidence that the RI-CLPM, just like the DPM, are not equipped
# to handle time-varying effects of baseline confounders. In this script we do the following:
# 1) The user to passes in B-matrices directly:
#       - Either a single 2 x k matrix (recycled across T waves), OR
#       - A list of length T containing 2 x k matrices (one per wave), OR
#       - A named list of scenarios, where each scenario contains either of the above.
# 2) Instead of generating B trajectories internally (constant/stepwise), the engine uses the user-supplied B lists.
# 3) We simulate panel data under a linear CLPM data-generating process with time-varying B matrices.
# 4) Build the models dynamically for adaption to T: CLPM, RI-CLPM, DPM, CLPM with confounders, and linear BCA in conjunction with CLPM. 
# 5) For all extracted lagged parameters (AR and cross-lags), we now also return p-values and confidence intervals.

library(mvtnorm)
library(lavaan)
library(tidyverse)   
library(parallel)
library(pbapply)

############################################################
##  1. User-supplied B-matrices helpers 
############################################################

# helper: coerce a user-provided B input into a list of length T
# accepts either:
#   - a single 2 x k matrix (recycled to all waves), or
#   - a list of length T, each element a 2 x k matrix
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
      stop("B list must have length T (one 2 x k matrix per wave).")

    for (t in 1:T) {
      if (!is.matrix(B_in[[t]]))
        stop("Each element of B_list must be a 2 x k matrix.")
      if (nrow(B_in[[t]]) != 2)
        stop("Each B_t must have exactly 2 rows (X, Y).")
    }

    # set names if missing
    if (is.null(names(B_in))) names(B_in) <- paste0("t", 1:T)

    return(B_in)
  }

  stop("B_in must be either a 2 x k matrix or a list of length T of 2 x k matrices.")
}

############################################################
##  3. Simulate panel data
############################################################

# this function simulates panel data under a CLPM data-generating process, where 
# 1) we simulate the confounders first for each wave, based on the given Psi covariance matrix
# 2) we compute the variance these confounders induce at each wave given the B matrix for that wave
# 3) we compute how much variance must then be induced by the dynamic process to reach the target covariance 
# 4) however, we want some extra covariance between X and Y at each occasion, meaning that we need to first compute
#    the implied covariance at each wave, and then add this extra covariance to the target covariance at each wave
# 5) calculate the innovations covariance needed to reach this target covariance at each wave
# 6) simulate the panel data panel data for each wave using these innovations
#    then adding the lagged effects from A and the direct confounder effects from B
# 7) repeat for T waves with varying B matrices

simulate_panel_data <- function(
    N,                                                         # number of individuals
    T,                                                         # number of waves
    A,                                                         # 2x2 autoregressive/cross-lag matrix
    B_list,                                                    # list of B matrices: B_list[[t]] is 2 x k
    Psi,                                                       # k x k confounder covariance matrix
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

  # number of confounders
  k <- ncol(Psi)

  # simulate the confounders from their multivariate normal distribution
  # given a specified covariance matrix Psi, where diag=1 and for now the off-diag=0, giving I, the identity matrix. 
  U <- mvtnorm::rmvnorm(
    n     = N,                                                 # number of individuals
    mean  = rep(0, k),                                         # all confounders have mean 0
    sigma = Psi                                                # covariance matrix of confounders
  )

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
    S_U_t <- B_t %*% Psi %*% t(B_t)

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
  df <- matrix(NA, nrow = N, ncol = 2*T + k)

  # set column names
  colnames(df) <- c(paste0("x", 1:T),
                    paste0("y", 1:T),
                    paste0("c", 1:k))

  # add confounders to dataframe
  df[, (2*T + 1):(2*T + k)] <- U

  # simulate the first wave, which is different because there are no lagged values yet
  Ddyn <- mvtnorm::rmvnorm(
    n     = N,                                                   # number of individuals
    mean  = c(0, 0),                                             # mean 0 for X and Y
    sigma = S_dyn_list[[1]]                                      # variance covariance matrix at wave 1
  ) 

  # add the direct confounder effects
  obs1 <- Ddyn + U %*% t(B_list[[1]])

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
    obs <- Ddyn + U %*% t(B_list[[t]])

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
##  5. Residualizer (linear confounders)
############################################################

residualise_panel_linearC <- function(df,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  
  # convert to data frame
  df <- as.data.frame(df)
  
  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value=TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value=TRUE)
  c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value=TRUE)

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

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

# same as above but for CLPM with residualised confounders
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

############################################################
##  7. Extract lagged parameters + rho (NOW incl. p-values + CI)
############################################################

# since all fits are slightly different, we need to figure out where in the model object
# our parameters of interest: autoregressive, cross-lagged, and residual correlation (rho) are located

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

############################################################
##  8. Wrapper for one replication
############################################################

# so now we need to obtain a function that runs all other functions for one replication of the study
# this function should contain all the arguments that are defined in the previous sections

# Instead of sampling/generating B, the user passes B_scenarios:
#   - named list: B_scenarios[["constant"]] = (matrix or list length T)
#   - or unnamed list aligned with `scenarios` order
run_one_rep_study1 <- function(
    rep_id,                                                 # replication index (set by outer loop)
    N,                                                      # sample size
    T,                                                      # number of waves
    k,                                                      # number of confounders
    scenarios,                                              # character vector: e.g. c("constant","stepwise",...)
    B_scenarios,                                            # user-supplied B matrices per scenario
    A,                                                      # 2×2 autoregressive + cross-lag matrix
    Psi,                                                    # k×k confounder covariance
    rho_extra,                                              # extra covariance added to X,Y each wave
    models_to_run,                                          # e.g. c("clpm","riclpm","dpm","lbca","adj")
    base_seed = 1234,                                       # base seed
    ci_level = 0.95                                         # CI level for extracted parameters
){
  
  # set seed for this replication, use rep_id so it varies
  set.seed(base_seed + rep_id)

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
    if (ncol(B_list[[1]]) != k)
      stop("Mismatch: user-supplied B has k = ", ncol(B_list[[1]]), " columns, but function argument k = ", k, ".")
    if (!all(dim(Psi) == c(k, k)))
      stop("Mismatch: Psi must be k x k with k = ", k, ".")

    # extract mean betas
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec

    # try to simulate panel data
    df <- tryCatch(
      simulate_panel_data(
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

        # cross-lag YX (est + p + CI)
        estYX_CLPM      = NA, pYX_CLPM      = NA, ciL_YX_CLPM      = NA, ciU_YX_CLPM      = NA,
        estYX_RI_CLPM   = NA, pYX_RI_CLPM   = NA, ciL_YX_RI_CLPM   = NA, ciU_YX_RI_CLPM   = NA,
        estYX_DPM       = NA, pYX_DPM       = NA, ciL_YX_DPM       = NA, ciU_YX_DPM       = NA,
        estYX_CLPM_Adj  = NA, pYX_CLPM_Adj  = NA, ciL_YX_CLPM_Adj  = NA, ciU_YX_CLPM_Adj  = NA,
        estYX_CLPM_LBCA = NA, pYX_CLPM_LBCA = NA, ciL_YX_CLPM_LBCA = NA, ciU_YX_CLPM_LBCA = NA,

        # AR X (est + p + CI)
        estA_CLPM      = NA, pA_CLPM      = NA, ciL_A_CLPM      = NA, ciU_A_CLPM      = NA,
        estA_RI_CLPM   = NA, pA_RI_CLPM   = NA, ciL_A_RI_CLPM   = NA, ciU_A_RI_CLPM   = NA,
        estA_DPM       = NA, pA_DPM       = NA, ciL_A_DPM       = NA, ciU_A_DPM       = NA,
        estA_CLPM_Adj  = NA, pA_CLPM_Adj  = NA, ciL_A_CLPM_Adj  = NA, ciU_A_CLPM_Adj  = NA,
        estA_CLPM_LBCA = NA, pA_CLPM_LBCA = NA, ciL_A_CLPM_LBCA = NA, ciU_A_CLPM_LBCA = NA,

        # AR Y (est + p + CI)
        estAY_CLPM      = NA, pAY_CLPM      = NA, ciL_AY_CLPM      = NA, ciU_AY_CLPM      = NA,
        estAY_RI_CLPM   = NA, pAY_RI_CLPM   = NA, ciL_AY_RI_CLPM   = NA, ciU_AY_RI_CLPM   = NA,
        estAY_DPM       = NA, pAY_DPM       = NA, ciL_AY_DPM       = NA, ciU_AY_DPM       = NA,
        estAY_CLPM_Adj  = NA, pAY_CLPM_Adj  = NA, ciL_AY_CLPM_Adj  = NA, ciU_AY_CLPM_Adj  = NA,
        estAY_CLPM_LBCA = NA, pAY_CLPM_LBCA = NA, ciL_AY_CLPM_LBCA = NA, ciU_AY_CLPM_LBCA = NA,

        estRho_CLPM = NA,

        fail_CLPM      = TRUE,
        fail_RI_CLPM   = TRUE,
        fail_DPM       = TRUE,
        fail_CLPM_Adj  = TRUE,
        fail_CLPM_LBCA = TRUE,

        err_CLPM      = "sim failed",
        err_RI_CLPM   = "sim failed",
        err_DPM       = "sim failed",
        err_CLPM_Adj  = "sim failed",
        err_CLPM_LBCA = "sim failed",

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
    res_lbca <- if ("lbca"   %in% models_to_run) safe_fit_clpm_resid(model_clpm, df) else list(fit=NULL, err=NA)

    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit

    # extract lagged parameters (NOW includes p-values + CI)
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",   ci_level = ci_level)
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm", ci_level = ci_level)
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm",    ci_level = ci_level)
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm",   ci_level = ci_level)
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm",   ci_level = ci_level)

    # residual correlations
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm")

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

      # residual correlation
      estRho_CLPM = rho_clpm,

      # failure indicators
      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),

      # error messages
      err_CLPM      = rep(res_clpm$err,   T),
      err_RI_CLPM   = rep(res_ric$err,    T),
      err_DPM       = rep(res_dpm0$err,   T),
      err_CLPM_Adj  = rep(res_adj$err,    T),
      err_CLPM_LBCA = rep(res_lbca$err,   T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy$est, lag_ric$xy$est, lag_dpm0$xy$est, lag_adj$xy$est, lag_lbca$xy$est
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}

############################################################
## 9. Main simulation function — PARALLEL
############################################################

run_simulation_study1 <- function(
    reps,                                                                    # number of replications
    N,                                                                       # sample size
    T,                                                                       # number of waves
    k,                                                                       # number of confounders
    scenarios,                                                               # e.g., c("constant","stepwise")
    B_scenarios,                                                             # user-supplied B matrices per scenario
    A,                                                                       # 2×2 AR + cross-lag matrix
    Psi,                                                                     # k×k confounder covariance
    rho_extra,                                                               # extra covariance added to observations
    models_to_run,                                                           # c("clpm","riclpm","dpm","adj","lbca")
    cores = NULL,                                                            # default is detectCores()/2
    base_seed = 1234,                                                        # master seed for reproducible reps
    ci_level = 0.95                                                          # CI level for extracted parameters
) {

  # if the number of cores is not specified, detect and use half of available cores
  if (is.null(cores)) {

    # detect and use half of available cores
    cores <- max(1, floor(parallel::detectCores() / 2))
  }

  # if cores is 1, run sequentially without parallelization
  if (cores == 1L) {

    # run sequentially
    results_list <- lapply(
      X = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study1(
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
    NULL
  })

  # export all necessary functions and variables to the cluster
  parallel::clusterExport(
    cl,
    c(
      "as_B_list",
      "simulate_panel_data",
      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",
      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",
      "extract_lagged_parameters",
      "extract_rho_vec",
      "residualise_panel_linearC",
      "run_one_rep_study1",
      "N","T","k","scenarios",
      "B_scenarios",
      "A","Psi","rho_extra","models_to_run","base_seed","ci_level"
    ),
    envir = environment()
  )

  # run the simulation with a progress bar
  results_list <- pbapply::pblapply(
    X = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study1(
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
