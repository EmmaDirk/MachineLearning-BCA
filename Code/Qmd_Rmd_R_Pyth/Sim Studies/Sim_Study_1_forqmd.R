############################################################
##  SIMULATION STUDY 1 — COMPLETE FUNCTION SCRIPT
##  (PARALLEL, WITH RHO, LONG FORMAT OUTPUT)
############################################################

############################################################
##  LIBRARIES
############################################################

library(mvtnorm)
library(lavaan)
library(tidyverse)   # for bind_rows etc. (not used on workers)
library(parallel)
library(pbapply)

############################################################
##  1. Sample baseline B-matrix (linear confounders)
############################################################

sample_B_linear <- function(
    k,
    R2_1,
    min_abs   = 0.01,
    max_abs   = 0.30,
    max_tries = 100000
) {

  target_X <- R2_1 / 2
  target_Y <- R2_1 / 2

  for (i in seq_len(max_tries)) {

    u_x <- rnorm(k)
    u_x <- u_x / sqrt(sum(u_x^2))
    b_x <- sqrt(target_X) * u_x

    u_y <- rnorm(k)
    u_y <- u_y / sqrt(sum(u_y^2))
    b_y <- sqrt(target_Y) * u_y

    if (all(abs(b_x) >= min_abs,
            abs(b_x) <= max_abs,
            abs(b_y) >= min_abs,
            abs(b_y) <= max_abs)) {

      B1 <- rbind(b_x, b_y)
      rownames(B1) <- c("X", "Y")
      colnames(B1) <- paste0("c", 1:k)
      return(B1)
    }
  }

  stop("Failed to sample a valid B matrix within max_tries.")
}

############################################################
##  2. Generate B trajectory (same substantive logic)
############################################################

generate_B_trajectory <- function(
    B1,
    T,
    scenario = c("constant","linear","sinusoidal","stepwise","random_walk"),
    target_sd = 0.10,
    rw_sd     = 0.05
){

  scenario <- match.arg(scenario)
  k <- ncol(B1)

  if (scenario == "constant") {
    v <- rep(0, T)

  } else if (scenario == "linear") {
    v <- seq(0, 1, length.out = T)

  } else if (scenario == "sinusoidal") {
    v <- sin(seq(0, 2*pi, length.out = T))

  } else if (scenario == "stepwise") {
    v <- c(rep(0, floor(T/2)), rep(1, T - floor(T/2)))

  } else if (scenario == "random_walk") {
    steps <- rnorm(T - 1, mean = 0, sd = rw_sd)
    v <- c(0, cumsum(steps))
  }

  v_centered <- v - mean(v)
  current_sd <- sd(v_centered)

  if (current_sd > 0) {
    scaling_factor <- target_sd / current_sd
    v_scaled <- v_centered * scaling_factor
  } else {
    v_scaled <- rep(0, T)
  }

  B_list <- vector("list", T)
  names(B_list) <- paste0("t", 1:T)

  for (t in 1:T) {
    B_list[[t]] <- B1 + v_scaled[t]
  }

  B_list
}

############################################################
##  3. Simulate panel data
############################################################

simulate_panel_data <- function(
    N,
    T,
    B_list,
    ar    = 0.25,
    cross = 0.10,
    rho   = 0.30
){

  # 2x2 transition matrix
  A <- matrix(c(ar, cross,
                cross, ar), 2, 2, byrow = TRUE)

  k <- ncol(B_list[[1]])

  # linear confounders
  U <- mvtnorm::rmvnorm(N, sigma = diag(k))

  # container: x1..xT, y1..yT, c1..ck
  df <- matrix(NA, N, 2*T + k)
  colnames(df) <- c(paste0("x",1:T), paste0("y",1:T), paste0("c",1:k))
  df[, (2*T+1):(2*T+k)] <- U

  # helper to build innovation covariance so that total var ≈ 1
  make_Sigma_e <- function(Bt, S_dyn, rho){

    # confounder-induced covariance
    S_U <- Bt %*% t(Bt)

    # we want diag(S_dyn + S_U + Sigma_e) ≈ 1
    d <- 1 - diag(S_dyn + S_U)
    d[d < 1e-12] <- 1e-12

    R <- matrix(c(1, rho,
                  rho, 1), 2, 2)

    D <- diag(sqrt(d))
    Sigma_e <- D %*% R %*% D
    Sigma_e
  }

  # time 1
  S_dyn1   <- matrix(0, 2, 2)
  Sigma_e1 <- make_Sigma_e(B_list[[1]], S_dyn1, rho)

  Ddyn <- mvtnorm::rmvnorm(N, sigma = Sigma_e1)
  obs1 <- Ddyn + U %*% t(B_list[[1]])

  df[,1]   <- obs1[,1]
  df[,1+T] <- obs1[,2]

  S_prev <- Sigma_e1

  # later times
  for(t in 2:T){

    S_dyn_t   <- A %*% S_prev %*% t(A)
    Sigma_et  <- make_Sigma_e(B_list[[t]], S_dyn_t, rho)

    Ddyn <- Ddyn %*% t(A) + mvtnorm::rmvnorm(N, sigma = Sigma_et)

    obs <- Ddyn + U %*% t(B_list[[t]])
    df[, t]    <- obs[,1]
    df[, t+T ] <- obs[,2]

    S_prev <- S_dyn_t + Sigma_et
    S_prev <- (S_prev + t(S_prev))/2
  }

  df
}

############################################################
##  4. Model builders
############################################################

build_clpm <- function(T) {

  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

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

build_clpm_with_Cs <- function(T, k) {

  C_names <- paste0("c", 1:k, collapse=" + ")

  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(
        sprintf("x%d ~ x%d + y%d + %s", t, t-1, t-1, C_names),
        sprintf("y%d ~ x%d + y%d + %s", t, t-1, t-1, C_names)
      )
    })), collapse="\n"
  )

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

build_riclpm <- function(T) {

  ri_block <- paste0(
    "rix =~ ", paste(sprintf("1*x%d", 1:T), collapse=" + "), "\n",
    "riy =~ ", paste(sprintf("1*y%d", 1:T), collapse=" + "), "\n",
    "rix ~~ rix\n riy ~~ riy\n rix ~~ riy\n"
  )

  resid_fix <- paste0(
    paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse="; "), "\n",
    paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  within_lat <- paste0(
    paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse="; "), "\n",
    paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  orth <- paste0(
    "rix ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "rix ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n"
  )

  within_var <- paste0(
    paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n",
    paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse="; "), "\n"
  )

  within_cov <- paste0(
    paste(sprintf("wy%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n"
  )

  regress <- paste(
    unlist(lapply(2:T, function(t){
      c(
        sprintf("wx%d ~ wx%d + wy%d", t, t-1, t-1),
        sprintf("wy%d ~ wx%d + wy%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  means <- paste0(
    paste(paste0("x",1:T), collapse=" + "), " ~ mx*1\n",
    paste(paste0("y",1:T), collapse=" + "), " ~ my*1\n"
  )

  paste(ri_block, resid_fix, within_lat, orth,
        within_var, within_cov, regress, means, sep="\n")
}

build_dpm <- function(T) {

  FX_block <- paste0(
    "FX =~ ", paste(sprintf("1*x%d", 2:T), collapse=" + "), "\n"
  )

  FY_block <- paste0(
    "FY =~ ", paste(sprintf("1*y%d", 2:T), collapse=" + "), "\n"
  )

  fx_cov_block <- "FX ~~ x1 + y1\n"
  fy_cov_block <- "FY ~~ x1 + y1\n"

  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  resid_cov_block <- paste(
    sprintf("x%d ~~ y%d", 1:T, 1:T),
    collapse="\n"
  )

  latent_cov_block <- paste(
    "FX ~~ FX",
    "FY ~~ FY",
    "FX ~~ FY",
    sep="\n"
  )

  resid_var_block <- paste(
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  means_block <- paste(
    paste(sprintf("x%d", 1:T), collapse=" + "), "~ 1\n",
    paste(sprintf("y%d", 1:T), collapse=" + "), "~ 1\n"
  )

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

  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value=TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value=TRUE)
  c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value=TRUE)

  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  C <- as.matrix(df[c_cols])

  for (x in x_cols)
    df[[x]] <- resid(lm(df[[x]] ~ C))

  for (y in y_cols)
    df[[y]] <- resid(lm(df[[y]] ~ C))

  df
}

############################################################
##  6. SAFE FITTING HELPERS (capture error messages)
############################################################

safe_fit_clpm <- function(model_string, data) {

  err <- NA_character_

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "MLR",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(fit = fit, err = err)
}

safe_fit_riclpm <- function(model_string, data) {

  err <- NA_character_

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "MLR",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(fit = fit, err = err)
}

safe_fit_dpm <- function(model_string, data) {

  err <- NA_character_

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "MLR",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(fit = fit, err = err)
}

safe_fit_clpm_C <- function(model_string, data) {

  err <- NA_character_

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = as.data.frame(data),
      estimator = "MLR",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(fit = fit, err = err)
}

safe_fit_clpm_resid <- function(model_string, data) {

  err <- NA_character_

  df_resid <- tryCatch(
    residualise_panel_linearC(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_resid)) {
    return(list(fit = NULL, err = err))
  }

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_resid,
      estimator = "MLR",
      warn      = FALSE
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )

  list(fit = fit, err = err)
}

############################################################
##  7. Extract lagged parameters + rho
############################################################

extract_lagged_parameters <- function(
    fit,
    T,
    model_type = c("clpm", "riclpm", "dpm")
){
  model_type <- match.arg(model_type)

  if (is.null(fit)) {
    return(list(
      ar_x = rep(NA, T-1),
      ar_y = rep(NA, T-1),
      xy   = rep(NA, T-1),
      yx   = rep(NA, T-1)
    ))
  }

  pe <- tryCatch(lavaan::parameterEstimates(fit), error=function(e) NULL)
  if (is.null(pe)) {
    return(list(
      ar_x = rep(NA, T-1),
      ar_y = rep(NA, T-1),
      xy   = rep(NA, T-1),
      yx   = rep(NA, T-1)
    ))
  }

  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  grab <- function(lhs, rhs){
    ix <- which(pe$lhs == lhs & pe$rhs == rhs)
    if (length(ix) == 0) return(NA_real_)
    pe$est[ix[1]]
  }

  ar_x <- numeric(T-1)
  ar_y <- numeric(T-1)
  xy   <- numeric(T-1)
  yx   <- numeric(T-1)

  for (t in 2:T) {
    ar_x[t-1] <- grab(paste0(xvar, t), paste0(xvar, t-1))
    ar_y[t-1] <- grab(paste0(yvar, t), paste0(yvar, t-1))
    xy[t-1]   <- grab(paste0(yvar, t), paste0(xvar, t-1))
    yx[t-1]   <- grab(paste0(xvar, t), paste0(yvar, t-1))
  }

  list(
    ar_x = ar_x,
    ar_y = ar_y,
    xy   = xy,
    yx   = yx
  )
}

# Residual correlation (rho) per occasion from CLPM/others
extract_rho_vec <- function(
    fit,
    T,
    model_type = c("clpm","riclpm","dpm")
){
  model_type <- match.arg(model_type)

  if (is.null(fit)) return(rep(NA_real_, T))

  pe <- tryCatch(lavaan::parameterEstimates(fit), error=function(e) NULL)
  if (is.null(pe)) return(rep(NA_real_, T))

  if (model_type == "riclpm") {
    xvar <- "wx"
    yvar <- "wy"
  } else {
    xvar <- "x"
    yvar <- "y"
  }

  rho <- numeric(T)

  for (t in 1:T) {
    lhs_xy <- paste0(xvar, t)
    lhs_yx <- paste0(yvar, t)

    ix <- which(pe$lhs == lhs_xy & pe$rhs == lhs_yx)
    if (length(ix) == 0) {
      ix <- which(pe$lhs == lhs_yx & pe$rhs == lhs_xy)
    }

    if (length(ix) == 0) {
      rho[t] <- NA_real_
      next
    }

    cov_xy <- pe$est[ix[1]]

    vx_idx <- which(pe$lhs == lhs_xy & pe$rhs == lhs_xy)
    vy_idx <- which(pe$lhs == lhs_yx & pe$rhs == lhs_yx)

    if (length(vx_idx) == 0 || length(vy_idx) == 0) {
      rho[t] <- NA_real_
      next
    }

    vx <- pe$est[vx_idx[1]]
    vy <- pe$est[vy_idx[1]]

    if (is.na(vx) || is.na(vy) || vx <= 0 || vy <= 0) {
      rho[t] <- NA_real_
    } else {
      rho[t] <- cov_xy / sqrt(vx * vy)
    }
  }

  rho
}

############################################################
##  8. Wrapper for one replication
############################################################

run_one_rep_study1 <- function(
  rep_id,
  N,
  T,
  k,                 # number of linear confounders
  R2_1,              # total confounder R^2 at t = 1
  target_sd,         # SD of beta variation over time
  scenarios_internal,# character vector: "constant", "linear", ...
  scenarios_pretty,  # pretty labels for output
  ar,                # true AR (x_t <- x_{t-1})
  cross,             # true cross-lag (y_t <- x_{t-1})
  rho,               # innovation correlation
  model_clpm,
  model_riclpm,
  model_dpm,
  model_clpm_with_Cs,
  base_seed = 12345
){

  set.seed(base_seed + rep_id)

  # baseline B at t = 1
  B1 <- sample_B_linear(
    k    = k,
    R2_1 = R2_1
  )

  out_list <- vector("list", length(scenarios_internal))

  for (j in seq_along(scenarios_internal)) {

    scen_int    <- scenarios_internal[j]
    scen_pretty <- scenarios_pretty[j]

    # B trajectory for this scenario
    B_list <- generate_B_trajectory(
      B1        = B1,
      T         = T,
      scenario  = scen_int,
      target_sd = target_sd
    )

    # "beta" summaries (mean of row 1 and 2 over confounders)
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec  # same as X-row mean, as in your prints

    # simulate data
    df <- tryCatch(
      simulate_panel_data(
        N      = N,
        T      = T,
        B_list = B_list,
        ar     = ar,
        cross  = cross,
        rho    = rho
      ),
      error = function(e) NULL
    )

    # if simulation fails, return NA row for this scenario
    if (is.null(df)) {
      out_list[[j]] <- data.frame(
        run        = rep(rep_id, T),
        occasion   = 1:T,
        scenario   = scen_pretty,
        beta       = beta_vec,
        beta_X     = beta_X_vec,
        beta_Y     = beta_Y_vec,
        true_cross = cross,
        true_auto  = ar,

        est_CLPM        = NA_real_,
        est_RI_CLPM     = NA_real_,
        est_DPM         = NA_real_,
        est_CLPM_LBCA   = NA_real_,
        est_CLPM_Adj    = NA_real_,

        estA_CLPM       = NA_real_,
        estA_RI_CLPM    = NA_real_,
        estA_DPM        = NA_real_,
        estA_CLPM_LBCA  = NA_real_,
        estA_CLPM_Adj   = NA_real_,

        estRho_CLPM     = NA_real_,

        fail_CLPM       = TRUE,
        fail_RI_CLPM    = TRUE,
        fail_DPM        = TRUE,
        fail_CLPM_LBCA  = TRUE,
        fail_CLPM_Adj   = TRUE,

        err_CLPM        = "sim failed",
        err_RI_CLPM     = "sim failed",
        err_DPM         = "sim failed",
        err_CLPM_LBCA   = "sim failed",
        err_CLPM_Adj    = "sim failed",

        is_na_run       = 1L
      )
      next
    }

    ## fit models (safe) -----------------------------------

    res_clpm   <- safe_fit_clpm(model_clpm, df)
    res_ric    <- safe_fit_riclpm(model_riclpm, df)
    res_dpm    <- safe_fit_dpm(model_dpm, df)
    res_adj    <- safe_fit_clpm_C(model_clpm_with_Cs, df)
    res_lbca   <- safe_fit_clpm_resid(model_clpm, df)

    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit

    ## extract lagged parameters ---------------------------

    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm")
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm")
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm")
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm")
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm")

    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm")

    fail_CLPM      <- is.null(fit_clpm_raw)
    fail_RI_CLPM   <- is.null(fit_ric)
    fail_DPM       <- is.null(fit_dpm0)
    fail_CLPM_Adj  <- is.null(fit_adj)
    fail_CLPM_LBCA <- is.null(fit_lbca)

    # one error string per scenario, repeated over occasions
    err_CLPM      <- rep(res_clpm$err,   T)
    err_RI_CLPM   <- rep(res_ric$err,    T)
    err_DPM       <- rep(res_dpm$err,    T)
    err_CLPM_Adj  <- rep(res_adj$err,    T)
    err_CLPM_LBCA <- rep(res_lbca$err,   T)

    out_list[[j]] <- data.frame(
      run        = rep(rep_id, T),
      occasion   = 1:T,
      scenario   = scen_pretty,

      beta       = beta_vec,
      beta_X     = beta_X_vec,
      beta_Y     = beta_Y_vec,

      true_cross = cross,
      true_auto  = ar,

      est_CLPM      = c(NA, lag_raw$xy),
      est_RI_CLPM   = c(NA, lag_ric$xy),
      est_DPM       = c(NA, lag_dpm0$xy),
      est_CLPM_LBCA = c(NA, lag_lbca$xy),
      est_CLPM_Adj  = c(NA, lag_adj$xy),

      estA_CLPM      = c(NA, lag_raw$ar_x),
      estA_RI_CLPM   = c(NA, lag_ric$ar_x),
      estA_DPM       = c(NA, lag_dpm0$ar_x),
      estA_CLPM_LBCA = c(NA, lag_lbca$ar_x),
      estA_CLPM_Adj  = c(NA, lag_adj$ar_x),

      estRho_CLPM    = rho_clpm,

      fail_CLPM      = fail_CLPM,
      fail_RI_CLPM   = fail_RI_CLPM,
      fail_DPM       = fail_DPM,
      fail_CLPM_LBCA = fail_CLPM_LBCA,
      fail_CLPM_Adj  = fail_CLPM_Adj,

      err_CLPM       = err_CLPM,
      err_RI_CLPM    = err_RI_CLPM,
      err_DPM        = err_DPM,
      err_CLPM_LBCA  = err_CLPM_LBCA,
      err_CLPM_Adj   = err_CLPM_Adj,

      is_na_run = as.integer(
        all(is.na(c(lag_raw$xy,
                    lag_ric$xy,
                    lag_dpm0$xy,
                    lag_adj$xy,
                    lag_lbca$xy)))
      )
    )
  }

  dplyr::bind_rows(out_list)
}

############################################################
## 9. Main simulation function — PARALLEL
############################################################

run_simulation_study1 <- function(
  reps,
  N,
  T,
  k,
  R2_1,
  target_sd,
  scenarios = c("Constant", "Linear", "Sinusoidal", "Stepwise", "Random Walk"),
  ar    = 0.25,
  cross = 0.10,
  rho   = 0.30,
  base_seed = 12345
){

  # model syntax
  model_clpm   <- build_clpm(T)
  model_riclpm <- build_riclpm(T)
  model_dpm    <- build_dpm(T)
  model_clpm_C <- build_clpm_with_Cs(T, k)

  # pretty vs internal scenario labels
  scenarios_pretty   <- scenarios
  scenarios_internal <- tolower(scenarios)
  scenarios_internal[scenarios_internal == "random walk"] <- "random_walk"

  # cluster setup
  cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(cores)

  parallel::clusterExport(
    cl,
    c(
      "sample_B_linear",
      "generate_B_trajectory",
      "simulate_panel_data",

      "build_clpm",
      "build_riclpm",
      "build_dpm",
      "build_clpm_with_Cs",

      "residualise_panel_linearC",

      "safe_fit_clpm",
      "safe_fit_riclpm",
      "safe_fit_dpm",
      "safe_fit_clpm_C",
      "safe_fit_clpm_resid",

      "extract_lagged_parameters",
      "extract_rho_vec",

      "run_one_rep_study1"
    ),
    envir = environment()
  )

  parallel::clusterEvalQ(cl, {
    library(lavaan)
    NULL
  })

  results_list <- pbapply::pblapply(
    X  = 1:reps,
    cl = cl,
    FUN = function(rep_id) {
      run_one_rep_study1(
        rep_id            = rep_id,
        N                 = N,
        T                 = T,
        k                 = k,
        R2_1              = R2_1,
        target_sd         = target_sd,
        scenarios_internal= scenarios_internal,
        scenarios_pretty  = scenarios_pretty,
        ar                = ar,
        cross             = cross,
        rho               = rho,
        model_clpm        = model_clpm,
        model_riclpm      = model_riclpm,
        model_dpm         = model_dpm,
        model_clpm_with_Cs= model_clpm_C,
        base_seed         = base_seed
      )
    }
  )

  parallel::stopCluster(cl)

  out_long <- dplyr::bind_rows(results_list)
  out_long
}

############################################################
## 10. Example run
############################################################

#set.seed(1)
#out_long <- run_simulation_study1(
#   reps      = 10,
#   N         = 1000,
#   T         = 5,
#   k         = 3,
#   R2_1      = 0.25,
#   target_sd = 0.10,
#   ar        = 0.25,
#   cross     = 0.10,
#   rho       = 0.30
#)




