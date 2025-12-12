##  STUDY 2 — NONLINEAR CONFOUNDING + ML RESIDUALIZERS
##  Extension of Study 1
##  - Arbitrary number of linear confounders k_lin
##  - One nonlinear confounder per linear confounder
##  - Total confounder R^2 at t=1: R2_1
##  - Nonlinear share at t=1: eta_1
##  - 5 time-varying scenarios (same SD of B-variation): 
##      Constant, Linear, Sinusoidal, Stepwise, Random Walk
##  - No shrink parameter
##  - Flexible cores for parallel execution
############################################################

############################################################
##  LIBRARIES
############################################################

library(mvtnorm)
library(lavaan)
library(dplyr)
library(tidyr)
library(parallel)
library(pbapply)

## progress bar options
pbapply::pboptions(type = "timer")

## ML residualizers
library(xgboost)
library(glmnet)
library(ranger)
library(dbarts)

############################################################
##  1. SAMPLE BASELINE B MATRIX AT t = 1
##  (linear + nonlinear confounders)
############################################################
## k_lin   : number of linear confounders
## R2_1    : total R^2 due to all confounders at t=1
## eta_1   : proportion of that strength in nonlinear part
############################################################

sample_B_matrix <- function(
    k_lin,
    R2_1, eta_1,
    min_abs   = 0.01,
    max_abs   = 0.30,
    max_tries = 5e5
){
  k_nonlin <- k_lin
  k        <- k_lin + k_nonlin

  ## shares
  LX <- (1 - eta_1) * R2_1 / 2
  LY <- (1 - eta_1) * R2_1 / 2
  NX <- eta_1       * R2_1 / 2
  NY <- eta_1       * R2_1 / 2

  for (i in 1:max_tries) {
    uxL <- rnorm(k_lin);    uxL <- uxL / sqrt(sum(uxL^2))
    uyL <- rnorm(k_lin);    uyL <- uyL / sqrt(sum(uyL^2))
    uxN <- rnorm(k_nonlin); uxN <- uxN / sqrt(sum(uxN^2))
    uyN <- rnorm(k_nonlin); uyN <- uyN / sqrt(sum(uyN^2))

    gamma_x <- c(sqrt(LX) * uxL, sqrt(NX) * uxN)
    gamma_y <- c(sqrt(LY) * uyL, sqrt(NY) * uyN)

    if (all(abs(gamma_x) >= min_abs & abs(gamma_x) <= max_abs &
            abs(gamma_y) >= min_abs & abs(gamma_y) <= max_abs)) {
      B0 <- rbind(gamma_x, gamma_y)
      rownames(B0) <- c("X", "Y")
      colnames(B0) <- c(
        paste0("c", 1:k_lin),
        paste0("c", 1:k_lin, "_NL")
      )
      return(B0)
    }
  }

  stop("Failed to sample B within max_tries.")
}

############################################################
##  2. TRAJECTORY GENERATOR — STUDY 1 LOGIC
############################################################
## target_sd = SD of B-variation across ALL scenarios
############################################################

generate_B_trajectory <- function(
    B1,
    T,
    scenario  = c("constant","linear","sinusoidal","stepwise","random_walk"),
    target_sd = 0.10,
    rw_sd     = 0.05
){
  scenario <- match.arg(scenario)

  if (scenario == "constant") {
    v <- rep(0, T)
  } else if (scenario == "linear") {
    v <- seq(0, 1, length.out = T)
  } else if (scenario == "sinusoidal") {
    v <- sin(seq(0, 2*pi, length.out = T))
  } else if (scenario == "stepwise") {
    v <- c(rep(0, floor(T/2)), rep(1, T - floor(T/2)))
  } else { # random_walk
    steps <- rnorm(T-1, mean = 0, sd = rw_sd)
    v <- c(0, cumsum(steps))
  }

  ## Center and rescale to common SD = target_sd
  v_centered <- v - mean(v)
  if (sd(v_centered) > 0) {
    v_scaled <- v_centered * target_sd / sd(v_centered)
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
##  3. ENFORCE SAME R^2 SHARE AT EVERY OCCASION
############################################################
## Ensures total confounder R^2 = R2_1 and nonlinear share = eta_1
## at every time point.
############################################################

enforce_R2_share <- function(
    B_list,
    R2_1, eta_1
){
  lapply(B_list, function(Bt) {
    lin_idx    <- grep("^c[0-9]+$", colnames(Bt))
    nonlin_idx <- grep("_NL$",    colnames(Bt))

    lin_norm    <- sum(Bt[1, lin_idx]^2)    + sum(Bt[2, lin_idx]^2)
    nonlin_norm <- sum(Bt[1, nonlin_idx]^2) + sum(Bt[2, nonlin_idx]^2)

    target_lin    <- (1 - eta_1) * R2_1
    target_nonlin <- eta_1       * R2_1

    Bt[, lin_idx]    <- Bt[, lin_idx] *
      sqrt(target_lin    / (lin_norm    + 1e-12))
    Bt[, nonlin_idx] <- Bt[, nonlin_idx] *
      sqrt(target_nonlin / (nonlin_norm + 1e-12))

    Bt
  })
}

############################################################
##  4. NONLINEAR PANEL SIMULATOR (USES B_list OVERRIDE)
############################################################

simulate_panel_nonlinear <- function(
    N, T,
    k_lin,
    ax, ay, bx, by,    # AR + cross-lag parameters
    rho,               # innovation correlation
    R2_1, eta_1,       # included for clarity/bookkeeping
    B_list_override,
    nonlin_type = c("cubic","exp","sin","tanh")
){

  nonlin_type <- match.arg(nonlin_type)

  ## 4.1 Linear confounders
  U_lin <- mvtnorm::rmvnorm(N, sigma = diag(k_lin))
  colnames(U_lin) <- paste0("c", 1:k_lin)

  ## 4.2 Nonlinear transforms (one per linear C)
  build_nonlin <- function(col){
    if (nonlin_type == "cubic") {
      return(scale(col^3))
    } else if (nonlin_type == "exp") {
      z <- exp(0.7 * col)
      return(scale(resid(lm(z ~ col))))
    } else if (nonlin_type == "sin") {
      z <- sin(3 * col)
      return(scale(resid(lm(z ~ col))))
    } else { # tanh
      z <- tanh(2 * col)
      return(scale(resid(lm(z ~ col))))
    }
  }

  U_nonlin <- do.call(
    cbind,
    lapply(1:k_lin, function(j) build_nonlin(U_lin[, j]))
  )
  colnames(U_nonlin) <- paste0("c", 1:k_lin, "_NL")

  U <- cbind(U_lin, U_nonlin)
  Psi <- cov(U)

  ## 4.3 B trajectory override
  Beta_list <- B_list_override
  k_total   <- ncol(U)

  ## 4.4 Dynamic matrix
  A <- matrix(c(ax, bx,
                by, ay), 2, 2, byrow = TRUE)

  ## helper: innovation covariance to keep Var(X_t), Var(Y_t) ≈ 1
  make_Sigma_e <- function(Bt, S_dyn, rho){
    S_U <- Bt %*% Psi %*% t(Bt)

    d <- 1 - diag(S_dyn + S_U)
    d[d < 1e-12] <- 1e-12

    R <- matrix(c(1, rho,
                  rho, 1), 2, 2)
    D <- diag(sqrt(d))
    D %*% R %*% D
  }

  ## 4.5 container
  df <- matrix(NA, N, 2*T + k_total)
  colnames(df) <- c(paste0("x",1:T), paste0("y",1:T), colnames(U))
  df[,(2*T+1):(2*T + k_total)] <- U

  ## 4.6 t = 1
  S_dyn   <- matrix(0, 2, 2)
  Sigma_e <- make_Sigma_e(Beta_list[[1]], S_dyn, rho)

  Ddyn <- mvtnorm::rmvnorm(N, sigma = Sigma_e)
  obs1 <- Ddyn + U %*% t(Beta_list[[1]])

  df[, 1]   <- obs1[,1]
  df[, 1+T] <- obs1[,2]

  S_prev <- Sigma_e

  ## 4.7 t >= 2
  for (t in 2:T) {
    S_dyn   <- A %*% S_prev %*% t(A)
    Sigma_e <- make_Sigma_e(Beta_list[[t]], S_dyn, rho)

    eps  <- mvtnorm::rmvnorm(N, sigma = Sigma_e)
    Ddyn <- Ddyn %*% t(A) + eps

    obs <- Ddyn + U %*% t(Beta_list[[t]])

    df[, t]    <- obs[, 1]
    df[, t+T ] <- obs[, 2]

    S_prev <- S_dyn + Sigma_e
    S_prev <- (S_prev + t(S_prev)) / 2
  }

  ## Summaries of B trajectory
  B_array <- tryCatch(simplify2array(Beta_list), error=function(e) NULL)
  if (!is.null(B_array)) {
    sd_B   <- mean(apply(B_array, c(1,2), sd, na.rm = TRUE), na.rm = TRUE)
    mean_B <- mean(B_array, na.rm = TRUE)
  } else {
    sd_B   <- NA_real_
    mean_B <- NA_real_
  }

  beta_X_vec <- sapply(Beta_list, function(Bt) mean(Bt[1, ], na.rm = TRUE))
  beta_Y_vec <- sapply(Beta_list, function(Bt) mean(Bt[2, ], na.rm = TRUE))

  list(
    df      = as.data.frame(df),
    B_list  = Beta_list,
    beta    = beta_X_vec,
    beta_X  = beta_X_vec,
    beta_Y  = beta_Y_vec,
    sd_B    = sd_B,
    mean_B  = mean_B
  )
}

############################################################
##  5. RESIDUALISERS
############################################################

residualise_panel_linearC <- function(df,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)

  ## Only pure linear Cs: c + digits only
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)
  if (length(c_cols) == 0) stop("No linear confounders.")

  for (x in x_cols) {
    df[[x]] <- resid(lm(df[[x]] ~ ., data = df[c_cols]))
  }
  for (y in y_cols) {
    df[[y]] <- resid(lm(df[[y]] ~ ., data = df[c_cols]))
  }

  df
}

residualise_panel_xgb <- function(df,
                                  x_prefix = "x",
                                  y_prefix = "y",
                                  c_prefix = "c",
                                  nrounds = 20,
                                  max_depth = 2,
                                  eta = 0.20) {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if (length(c_cols) == 0) stop("No linear confounders.")

  Xc <- as.matrix(df[c_cols])

  params <- list(
    objective = "reg:squarederror",
    eta       = eta,
    max_depth = max_depth,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    nthread   = 1
  )

  do_xgb <- function(y) {
    dtrain <- xgb.DMatrix(Xc, label = y)
    bst <- tryCatch(xgb.train(params, dtrain, nrounds, verbose = 0),
                    error = function(e) NULL)
    if (is.null(bst)) return(y)
    y - predict(bst, dtrain)
  }

  for (x in x_cols) df[[x]] <- do_xgb(df[[x]])
  for (y in y_cols) df[[y]] <- do_xgb(df[[y]])

  df
}

residualise_panel_lasso <- function(df,
                                    x_prefix = "x",
                                    y_prefix = "y",
                                    c_prefix = "c",
                                    alpha = 1) {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if (length(c_cols) == 0) stop("No linear confounders.")

  C <- scale(as.matrix(df[c_cols]))

  do_lasso <- function(y) {
    if (all(is.na(y))) return(y)
    cvfit <- tryCatch(
      glmnet::cv.glmnet(C, y, alpha = alpha),
      error = function(e) NULL
    )
    if (is.null(cvfit)) return(y)
    y_hat <- as.numeric(predict(cvfit, newx = C, s = "lambda.min"))
    y - y_hat
  }

  for (v in x_cols) df[[v]] <- do_lasso(df[[v]])
  for (v in y_cols) df[[v]] <- do_lasso(df[[v]])

  df
}

residualise_panel_ranger <- function(df,
                                     x_prefix = "x",
                                     y_prefix = "y",
                                     c_prefix = "c",
                                     num_trees = 200) {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if (length(c_cols) == 0) stop("No linear confounders.")

  Xc <- as.data.frame(df[c_cols])

  do_rf <- function(y) {
    dat <- cbind(Xc, y = y)
    fit <- tryCatch(
      ranger::ranger(
        y ~ .,
        data = dat,
        num.trees   = num_trees,
        write.forest = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(y)
    y_hat <- predict(fit, data = Xc)$predictions
    y - y_hat
  }

  for (v in x_cols) df[[v]] <- do_rf(df[[v]])
  for (v in y_cols) df[[v]] <- do_rf(df[[v]])

  df
}

residualise_panel_bart <- function(df,
                                   x_prefix = "x",
                                   y_prefix = "y",
                                   c_prefix = "c") {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if (length(c_cols) == 0) stop("No linear confounders.")

  Xc <- as.matrix(df[c_cols])

  do_bart <- function(y) {
    fit <- tryCatch(
      dbarts::bart(
        x.train   = Xc,
        y.train   = y,
        keeptrees = FALSE,
        verbose   = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$yhat.train.mean)) return(y)
    y_hat <- as.numeric(fit$yhat.train.mean)
    y - y_hat
  }

  for (v in x_cols) df[[v]] <- do_bart(df[[v]])
  for (v in y_cols) df[[v]] <- do_bart(df[[v]])

  df
}

############################################################
##  6. MODEL BUILDERS (CLPM, RI-CLPM, DPM, CLPM+Cs)
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

build_clpm_with_Cs <- function(T, k_lin) {

  C_names <- paste0("c", 1:k_lin, collapse=" + ")

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
##  7. SAFE FITTERS (CLPM, RI-CLPM, DPM, + residualizers)
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

## LBCA: residualize on linear Cs, then CLPM
safe_fit_clpm_lbca <- function(model_string, data) {
  err <- NA_character_
  df_res <- tryCatch(
    residualise_panel_linearC(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_res)) return(list(fit = NULL, err = err))

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_res,
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

## Adjusted CLPM: include linear Cs as covariates in the CLPM
safe_fit_clpm_adj <- function(model_string, data) {
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

safe_fit_clpm_xgb <- function(model_string, data) {
  err <- NA_character_
  df_res <- tryCatch(
    residualise_panel_xgb(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_res)) return(list(fit = NULL, err = err))

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_res,
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

safe_fit_clpm_lasso <- function(model_string, data) {
  err <- NA_character_
  df_res <- tryCatch(
    residualise_panel_lasso(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_res)) return(list(fit = NULL, err = err))

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_res,
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

safe_fit_clpm_rf <- function(model_string, data) {
  err <- NA_character_
  df_res <- tryCatch(
    residualise_panel_ranger(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_res)) return(list(fit = NULL, err = err))

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_res,
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

safe_fit_clpm_bart <- function(model_string, data) {
  err <- NA_character_
  df_res <- tryCatch(
    residualise_panel_bart(data),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(df_res)) return(list(fit = NULL, err = err))

  fit <- tryCatch(
    lavaan::lavaan(
      model_string,
      data      = df_res,
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
##  8. EXTRACT LAGGED PARAMETERS + rho
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
##  9. ONE REPLICATION
############################################################

run_one_rep_study2 <- function(
  rep_id,
  N, T,
  k_lin,
  ax, ay, bx, by, rho,
  R2_1, eta_1,
  target_sd,
  rw_sd,
  scenarios_internal,
  scenarios_pretty,
  model_clpm,
  model_riclpm,
  model_dpm,
  model_clpm_C,
  base_seed = 12345
){

  set.seed(base_seed + rep_id)

  ## baseline B at t = 1 (nonlinear)
  B0 <- sample_B_matrix(
    k_lin = k_lin,
    R2_1  = R2_1,
    eta_1 = eta_1
  )

  ## one nonlinear type per replication
  nl_type <- sample(c("cubic","exp","sin","tanh"), 1)

  out_list <- vector("list", length(scenarios_internal))

  for (j in seq_along(scenarios_internal)) {

    scen_int    <- scenarios_internal[j]
    scen_pretty <- scenarios_pretty[j]

    ## trajectory for this scenario
    B_list_raw <- generate_B_trajectory(
      B1        = B0,
      T         = T,
      scenario  = scen_int,
      target_sd = target_sd,
      rw_sd     = rw_sd
    )

    ## enforce nonlinear share
    B_list <- enforce_R2_share(
      B_list = B_list_raw,
      R2_1   = R2_1,
      eta_1  = eta_1
    )

    ## simulate data
    sim <- simulate_panel_nonlinear(
      N               = N,
      T               = T,
      k_lin           = k_lin,
      ax              = ax,
      ay              = ay,
      bx              = bx,
      by              = by,
      rho             = rho,
      R2_1            = R2_1,
      eta_1           = eta_1,
      B_list_override = B_list,
      nonlin_type     = nl_type
    )

    df        <- sim$df
    beta_vec  <- sim$beta
    beta_X    <- sim$beta_X
    beta_Y    <- sim$beta_Y
    sd_B      <- sim$sd_B
    mean_B    <- sim$mean_B

    ## fit models
    res_clpm      <- safe_fit_clpm(model_clpm,   df)
    res_ric       <- safe_fit_riclpm(model_riclpm, df)
    res_dpm       <- safe_fit_dpm(model_dpm,     df)
    res_lbca      <- safe_fit_clpm_lbca(model_clpm, df)
    res_adj       <- safe_fit_clpm_adj(model_clpm_C, df)
    res_lasso     <- safe_fit_clpm_lasso(model_clpm, df)
    res_rf        <- safe_fit_clpm_rf(model_clpm, df)
    res_bart      <- safe_fit_clpm_bart(model_clpm, df)
    res_xgb       <- safe_fit_clpm_xgb(model_clpm, df)

    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm$fit
    fit_lbca     <- res_lbca$fit
    fit_adj      <- res_adj$fit
    fit_lasso    <- res_lasso$fit
    fit_rf       <- res_rf$fit
    fit_bart     <- res_bart$fit
    fit_xgb      <- res_xgb$fit

    lag_raw    <- extract_lagged_parameters(fit_clpm_raw, T, "clpm")
    lag_ric    <- extract_lagged_parameters(fit_ric,       T, "riclpm")
    lag_dpm0   <- extract_lagged_parameters(fit_dpm0,      T, "dpm")
    lag_lbca   <- extract_lagged_parameters(fit_lbca,      T, "clpm")
    lag_adj    <- extract_lagged_parameters(fit_adj,       T, "clpm")
    lag_lasso  <- extract_lagged_parameters(fit_lasso,     T, "clpm")
    lag_rf     <- extract_lagged_parameters(fit_rf,        T, "clpm")
    lag_bart   <- extract_lagged_parameters(fit_bart,      T, "clpm")
    lag_xgb    <- extract_lagged_parameters(fit_xgb,       T, "clpm")

    rho_clpm   <- extract_rho_vec(fit_clpm_raw, T, "clpm")

    fail_CLPM        <- is.null(fit_clpm_raw)
    fail_RI_CLPM     <- is.null(fit_ric)
    fail_DPM         <- is.null(fit_dpm0)
    fail_CLPM_LBCA   <- is.null(fit_lbca)
    fail_CLPM_Adj    <- is.null(fit_adj)
    fail_CLPM_LASSO  <- is.null(fit_lasso)
    fail_CLPM_RF     <- is.null(fit_rf)
    fail_CLPM_BART   <- is.null(fit_bart)
    fail_CLPM_XGB    <- is.null(fit_xgb)

    err_CLPM        <- rep(res_clpm$err,      T)
    err_RI_CLPM     <- rep(res_ric$err,       T)
    err_DPM         <- rep(res_dpm$err,       T)
    err_CLPM_LBCA   <- rep(res_lbca$err,      T)
    err_CLPM_Adj    <- rep(res_adj$err,       T)
    err_CLPM_LASSO  <- rep(res_lasso$err,     T)
    err_CLPM_RF     <- rep(res_rf$err,        T)
    err_CLPM_BART   <- rep(res_bart$err,      T)
    err_CLPM_XGB    <- rep(res_xgb$err,       T)

    all_xy <- c(
      lag_raw$xy,
      lag_ric$xy,
      lag_dpm0$xy,
      lag_lbca$xy,
      lag_adj$xy,
      lag_lasso$xy,
      lag_rf$xy,
      lag_bart$xy,
      lag_xgb$xy
    )

    out_list[[j]] <- data.frame(
      run        = rep(rep_id, T),
      occasion   = 1:T,
      scenario   = scen_pretty,
      nl_type    = nl_type,

      beta       = beta_vec,
      beta_X     = beta_X,
      beta_Y     = beta_Y,

      true_cross = by,
      true_auto  = ax,

      est_CLPM        = c(NA, lag_raw$xy),
      est_RI_CLPM     = c(NA, lag_ric$xy),
      est_DPM         = c(NA, lag_dpm0$xy),
      est_CLPM_LBCA   = c(NA, lag_lbca$xy),
      est_CLPM_Adj    = c(NA, lag_adj$xy),
      est_CLPM_LASSO  = c(NA, lag_lasso$xy),
      est_CLPM_RF     = c(NA, lag_rf$xy),
      est_CLPM_BART   = c(NA, lag_bart$xy),
      est_CLPM_XGB    = c(NA, lag_xgb$xy),

      estA_CLPM        = c(NA, lag_raw$ar_x),
      estA_RI_CLPM     = c(NA, lag_ric$ar_x),
      estA_DPM         = c(NA, lag_dpm0$ar_x),
      estA_CLPM_LBCA   = c(NA, lag_lbca$ar_x),
      estA_CLPM_Adj    = c(NA, lag_adj$ar_x),
      estA_CLPM_LASSO  = c(NA, lag_lasso$ar_x),
      estA_CLPM_RF     = c(NA, lag_rf$ar_x),
      estA_CLPM_BART   = c(NA, lag_bart$ar_x),
      estA_CLPM_XGB    = c(NA, lag_xgb$ar_x),

      estRho_CLPM      = rho_clpm,

      sd_B             = sd_B,
      mean_B           = mean_B,

      fail_CLPM        = fail_CLPM,
      fail_RI_CLPM     = fail_RI_CLPM,
      fail_DPM         = fail_DPM,
      fail_CLPM_LBCA   = fail_CLPM_LBCA,
      fail_CLPM_Adj    = fail_CLPM_Adj,
      fail_CLPM_LASSO  = fail_CLPM_LASSO,
      fail_CLPM_RF     = fail_CLPM_RF,
      fail_CLPM_BART   = fail_CLPM_BART,
      fail_CLPM_XGB    = fail_CLPM_XGB,

      err_CLPM         = err_CLPM,
      err_RI_CLPM      = err_RI_CLPM,
      err_DPM          = err_DPM,
      err_CLPM_LBCA    = err_CLPM_LBCA,
      err_CLPM_Adj     = err_CLPM_Adj,
      err_CLPM_LASSO   = err_CLPM_LASSO,
      err_CLPM_RF      = err_CLPM_RF,
      err_CLPM_BART    = err_CLPM_BART,
      err_CLPM_XGB     = err_CLPM_XGB,

      is_na_run        = as.integer(all(is.na(all_xy))),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(out_list)
}

############################################################
##  10. MAIN SIMULATION FUNCTION — FLEXIBLE cores
############################################################

run_simulation_study2 <- function(
  reps,
  N,
  T,
  k_lin,
  ax,
  ay,
  bx,
  by,
  rho,
  R2_1,
  eta_1,
  target_sd = 0.10,
  rw_sd     = 0.10,
  scenarios = c("Constant", "Linear", "Sinusoidal", "Stepwise", "Random Walk"),
  cores     = max(1, parallel::detectCores() - 1),
  base_seed = 12345
){

  ## model syntax
  model_clpm   <- build_clpm(T)
  model_riclpm <- build_riclpm(T)
  model_dpm    <- build_dpm(T)
  model_clpm_C <- build_clpm_with_Cs(T, k_lin)

  scenarios_pretty   <- scenarios
  scenarios_internal <- tolower(scenarios)
  scenarios_internal[scenarios_internal == "random walk"] <- "random_walk"

  ## progress bar style
  pbapply::pboptions(type = "timer")

  if (cores > 1) {
    cl <- parallel::makeCluster(cores)

    parallel::clusterExport(
      cl,
      c(
        "sample_B_matrix",
        "generate_B_trajectory",
        "enforce_R2_share",
        "simulate_panel_nonlinear",

        "residualise_panel_linearC",
        "residualise_panel_xgb",
        "residualise_panel_lasso",
        "residualise_panel_ranger",
        "residualise_panel_bart",

        "build_clpm",
        "build_clpm_with_Cs",
        "build_riclpm",
        "build_dpm",

        "safe_fit_clpm",
        "safe_fit_riclpm",
        "safe_fit_dpm",
        "safe_fit_clpm_lbca",
        "safe_fit_clpm_adj",
        "safe_fit_clpm_xgb",
        "safe_fit_clpm_lasso",
        "safe_fit_clpm_rf",
        "safe_fit_clpm_bart",

        "extract_lagged_parameters",
        "extract_rho_vec",

        "run_one_rep_study2",

        "N","T","k_lin","ax","ay","bx","by","rho",
        "R2_1","eta_1","target_sd","rw_sd",
        "model_clpm","model_riclpm","model_dpm","model_clpm_C",
        "scenarios_internal","scenarios_pretty",
        "base_seed"
      ),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      library(mvtnorm)
      library(lavaan)
      library(xgboost)
      library(glmnet)
      library(ranger)
      library(dbarts)
      library(pbapply)
      NULL
    })

    results_list <- pbapply::pblapply(
      X  = 1:reps,
      cl = cl,
      FUN = function(rep_id) {
        run_one_rep_study2(
          rep_id            = rep_id,
          N                 = N,
          T                 = T,
          k_lin             = k_lin,
          ax                = ax,
          ay                = ay,
          bx                = bx,
          by                = by,
          rho               = rho,
          R2_1              = R2_1,
          eta_1             = eta_1,
          target_sd         = target_sd,
          rw_sd             = rw_sd,
          scenarios_internal= scenarios_internal,
          scenarios_pretty  = scenarios_pretty,
          model_clpm        = model_clpm,
          model_riclpm      = model_riclpm,
          model_dpm         = model_dpm,
          model_clpm_C      = model_clpm_C,
          base_seed         = base_seed
        )
      }
    )

    parallel::stopCluster(cl)

  } else {
    ## single-core run with progress bar
    results_list <- pbapply::pblapply(
      X  = 1:reps,
      FUN = function(rep_id) {
        run_one_rep_study2(
          rep_id            = rep_id,
          N                 = N,
          T                 = T,
          k_lin             = k_lin,
          ax                = ax,
          ay                = ay,
          bx                = bx,
          by                = by,
          rho               = rho,
          R2_1              = R2_1,
          eta_1             = eta_1,
          target_sd         = target_sd,
          rw_sd             = rw_sd,
          scenarios_internal= scenarios_internal,
          scenarios_pretty  = scenarios_pretty,
          model_clpm        = model_clpm,
          model_riclpm      = model_riclpm,
          model_dpm         = model_dpm,
          model_clpm_C      = model_clpm_C,
          base_seed         = base_seed
        )
      }
    )
  }

  out_long <- dplyr::bind_rows(results_list)
  out_long
}
