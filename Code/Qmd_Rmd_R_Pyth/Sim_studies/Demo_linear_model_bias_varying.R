############################################################
##  Libraries
############################################################

library(mvtnorm)
library(lavaan)
library(dplyr)
library(ggplot2)
library(parallel)
library(xgboost)
library(glmnet)

############################################################
##  1. Helper: sample B matrix given s0 and eta0
############################################################
# B is 2 x k: first row gamma_x, second row gamma_y
# k_lin    = number of linear confounders
# k_nonlin = number of nonlinear components (per confounder type)
# s0       = total confounder strength (approx Var_X^U + Var_Y^U)
# eta0     = fraction of that strength due to nonlinear components

sample_B_matrix <- function(
    k_lin,
    k_nonlin,
    s0,
    eta0,
    min_abs  = 0.01,
    max_abs  = 0.30,
    max_tries = 5e5
){
  k <- k_lin + k_nonlin

  # variance budgets for X and Y
  LX <- (1 - eta0) * s0 / 2
  LY <- (1 - eta0) * s0 / 2
  NX <- eta0 * s0 / 2
  NY <- eta0 * s0 / 2

  for (attempt in 1:max_tries) {

    # random directions (linear)
    uxL <- rnorm(k_lin)
    uyL <- rnorm(k_lin)
    # random directions (nonlinear)
    uxN <- rnorm(k_nonlin)
    uyN <- rnorm(k_nonlin)

    # normalize
    uxL <- uxL / sqrt(sum(uxL^2))
    uyL <- uyL / sqrt(sum(uyL^2))
    uxN <- uxN / sqrt(sum(uxN^2))
    uyN <- uyN / sqrt(sum(uyN^2))

    # scale to budgets
    gamma_x_lin <- sqrt(LX) * uxL
    gamma_y_lin <- sqrt(LY) * uyL
    gamma_x_non <- sqrt(NX) * uxN
    gamma_y_non <- sqrt(NY) * uyN

    gamma_x <- c(gamma_x_lin, gamma_x_non)
    gamma_y <- c(gamma_y_lin, gamma_y_non)

    if (
      all(abs(gamma_x) >= min_abs & abs(gamma_x) <= max_abs) &&
      all(abs(gamma_y) >= min_abs & abs(gamma_y) <= max_abs)
    ) {
      return(rbind(gamma_x, gamma_y))
    }
  }

  stop("No valid B matrix found within max_tries.")
}

############################################################
##  2. Residualisation functions
############################################################

# 2.1 Linear regression residualisation using ONLY linear confounders c1, c2, ...
residualise_panel_linearC <- function(df,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c"){
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if(length(c_cols) == 0){
    stop("No *linear* confounders (c1, c2, ...) found.")
  }

  for(x in x_cols){
    lm_fit <- lm(df[[x]] ~ ., data = df[c_cols])
    df[[x]] <- lm_fit$residuals
  }
  for(y in y_cols){
    lm_fit <- lm(df[[y]] ~ ., data = df[c_cols])
    df[[y]] <- lm_fit$residuals
  }

  df
}

# 2.2 XGBoost residualisation using ONLY linear confounders
residualise_panel_xgb <- function(df,
                                  x_prefix = "x",
                                  y_prefix = "y",
                                  c_prefix = "c",
                                  nrounds = 25,
                                  max_depth = 2,
                                  eta = 0.20,
                                  subsample = 0.8,
                                  colsample_bytree = 0.8){

  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_cols <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if(length(c_cols) == 0){
    stop("No *linear* confounders (c1, c2, ...) found.")
  }

  Xc <- as.matrix(df[, c_cols, drop = FALSE])

  params <- list(
    objective = "reg:squarederror",
    eta = eta,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    min_child_weight = 1,
    nthread = 1
  )

  for(x in x_cols){
    dtrain <- xgb.DMatrix(data = Xc, label = df[[x]])
    bst <- xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
    preds <- predict(bst, newdata = dtrain)
    df[[x]] <- df[[x]] - preds
  }

  for(y in y_cols){
    dtrain <- xgb.DMatrix(data = Xc, label = df[[y]])
    bst <- xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
    preds <- predict(bst, newdata = dtrain)
    df[[y]] <- df[[y]] - preds
  }

  df
}

# 2.3 LASSO residualisation: poly (up to 4) + pairwise interactions of linear C
# Uses only c1, c2, ... as input; nonlinear confounders are ignored.
residualise_panel_lasso <- function(df,
                                    x_prefix = "x",
                                    y_prefix = "y",
                                    c_prefix = "c",
                                    poly_order = 4,
                                    alpha = 1,
                                    nlambda = 200){

  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "[0-9]+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "[0-9]+$"), names(df), value = TRUE)
  c_linear <- grep(paste0("^", c_prefix, "[0-9]+$"), names(df), value = TRUE)

  if(length(c_linear) == 0){
    stop("No linear confounders detected (c1, c2, ...).")
  }

  C <- df[, c_linear, drop = FALSE]
  k_lin <- ncol(C)

  # polynomial terms
  poly_features <- lapply(seq_len(k_lin), function(j) {
    cj <- C[[j]]
    out <- lapply(1:poly_order, function(p) cj^p)
    do.call(cbind, out)
  })
  poly_mat <- do.call(cbind, poly_features)
  colnames(poly_mat) <- unlist(
    lapply(seq_len(k_lin), function(j){
      paste0(c_linear[j], "_p", 1:poly_order)
    })
  )

  # interactions
  inter_list <- list()
  inter_names <- c()
  idx <- 1
  if(k_lin > 1){
    for(i in 1:(k_lin-1)){
      for(j in (i+1):k_lin){
        inter_list[[idx]] <- C[[i]] * C[[j]]
        inter_names[idx] <- paste0(c_linear[i], "_x_", c_linear[j])
        idx <- idx + 1
      }
    }
  }

  if(length(inter_list) > 0){
    inter_mat <- do.call(cbind, inter_list)
    colnames(inter_mat) <- inter_names
  } else {
    inter_mat <- NULL
  }

  # full design matrix
  X_design <- cbind(C, poly_mat, inter_mat)
  X_design <- as.data.frame(X_design)

  # remove zero-variance columns
  nzv <- sapply(X_design, function(v) sd(v) > 1e-12)
  X_design <- X_design[, nzv, drop = FALSE]

  # convert to matrix
  X <- as.matrix(X_design)

  # helper: safe lasso fit
  safe_lasso <- function(y){
    fit <- tryCatch(
      glmnet(
        x = X,
        y = y,
        alpha = alpha,
        nlambda = nlambda,
        standardize = TRUE,
        intercept = TRUE
      ),
      error=function(e) NULL
    )

    if(is.null(fit)) return(y)  # fallback = no adjustment

    beta <- tryCatch(coef(fit, s = "lambda.min"), error=function(e) NULL)
    if(is.null(beta)) return(y)

    preds <- tryCatch(
      as.numeric(cbind(1, X) %*% beta),
      error=function(e) rep(0,length(y))
    )
    preds[is.na(preds)] <- 0

    y - preds
  }

  # residualise X variables
  for(x in x_cols){
    df[[x]] <- safe_lasso(df[[x]])
  }

  # residualise Y variables
  for(y in y_cols){
    df[[y]] <- safe_lasso(df[[y]])
  }

  df
}

############################################################
##  3. Model builders and fitters
############################################################

# 3.1 RI-CLPM (Hamaker style), using x1..xT, y1..yT
build_riclpm_model <- function(T){

  ri_block <- paste0(
    "rix =~ ", paste(paste0("1*x",1:T), collapse=" + "), "\n",
    "riy =~ ", paste(paste0("1*y",1:T), collapse=" + "), "\n",
    "rix ~~ rix\nriy ~~ riy\nrix ~~ riy\n"
  )

  resid_block <- paste0(
    paste(paste0("x",1:T," ~~ 0*x",1:T), collapse="; "), "\n",
    paste(paste0("y",1:T," ~~ 0*y",1:T), collapse="; "), "\n"
  )

  within_lat <- paste0(
    paste(paste0("wx",1:T," =~ 1*x",1:T), collapse="; "), "\n",
    paste(paste0("wy",1:T," =~ 1*y",1:T), collapse="; "), "\n"
  )

  orth <- paste0(
    "rix ~~ ", paste(paste0("0*wx",1:T), collapse=" + "), "\n",
    "rix ~~ ", paste(paste0("0*wy",1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(paste0("0*wx",1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(paste0("0*wy",1:T), collapse=" + "), "\n"
  )

  within_var <- paste0(
    paste(paste0("wx",1:T," ~~ wx",1:T), collapse="; "), "\n",
    paste(paste0("wy",1:T," ~~ wy",1:T), collapse="; "), "\n"
  )

  within_cov <- paste0(
    paste(paste0("wy",1:T," ~~ wx",1:T), collapse="; "), "\n"
  )

  regress <- paste(
    unlist(lapply(2:T, function(t){
      c(
        paste0("wx",t," ~ wx",t-1," + wy",t-1),
        paste0("wy",t," ~ wx",t-1," + wy",t-1)
      )
    })), collapse="\n"
  )

  means <- paste0(
    paste(paste0("x",1:T), collapse=" + "), " ~ mx*1\n",
    paste(paste0("y",1:T), collapse=" + "), " ~ my*1\n"
  )

  paste(
    ri_block,
    resid_block,
    within_lat,
    orth,
    within_var,
    within_cov,
    regress,
    means,
    sep="\n"
  )
}

fit_riclpm <- function(data, model_string){
  lavaan(model_string, data=as.data.frame(data), estimator="MLR")
}

# 3.2 CLPM (plain, no RI), using x1..xT, y1..yT
build_clpm_model <- function(T){

  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(
        paste0("x", t, " ~ x", t-1, " + y", t-1),
        paste0("y", t, " ~ x", t-1, " + y", t-1)
      )
    })), collapse="\n"
  )

  resid_cov <- paste(
    paste0("x", 1:T, " ~~ y", 1:T),
    collapse="\n"
  )

  resid_var <- paste(
    paste(paste0("x",1:T," ~~ x",1:T), collapse="\n"),
    paste(paste0("y",1:T," ~~ y",1:T), collapse="\n"),
    sep="\n"
  )

  means_block <- paste(
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  paste(
    regress_block, "\n",
    resid_cov, "\n",
    resid_var, "\n",
    means_block,
    sep="\n"
  )
}

fit_clpm <- function(data, model_string){
  lavaan(model_string, data = as.data.frame(data), estimator = "MLR")
}

# 3.3 DPM (dynamic panel model) with FX, FY
build_dpm_model <- function(T){

  FX_block <- paste0(
    "FX =~ ", paste(paste0("1*x", 2:T), collapse=" + "), "\n"
  )

  FY_block <- paste0(
    "FY =~ ", paste(paste0("1*y", 2:T), collapse=" + "), "\n"
  )

  FX_cov <- "FX ~~ x1 + y1\n"
  FY_cov <- "FY ~~ x1 + y1\n"

  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(
        paste0("x", t, " ~ x", t-1, " + y", t-1),
        paste0("y", t, " ~ x", t-1, " + y", t-1)
      )
    })), collapse="\n"
  )

  resid_cov <- paste(
    paste0("x", 1:T, " ~~ y", 1:T),
    collapse="\n"
  )

  latent_var <- paste(
    "FX ~~ FX",
    "FY ~~ FY",
    "FX ~~ FY",
    sep="\n"
  )

  resid_var <- paste(
    paste(paste0("x",1:T," ~~ x",1:T), collapse="\n"),
    paste(paste0("y",1:T," ~~ y",1:T), collapse="\n"),
    sep="\n"
  )

  means_block <- paste(
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  paste(
    FX_block,
    FY_block,
    FX_cov,
    FY_cov,
    regress_block, "\n",
    resid_cov, "\n",
    latent_var, "\n",
    resid_var, "\n",
    means_block,
    sep="\n"
  )
}

fit_dpm <- function(data, model_string){
  lavaan(model_string, data = as.data.frame(data), estimator = "MLR")
}

############################################################
##  4. Extract mean X->Y cross-lag for each model type
############################################################

extract_crosslag_mean <- function(fit, T, model_type = c("clpm","dpm","riclpm")){
  model_type <- match.arg(model_type)
  pe <- parameterEstimates(fit)

  if(model_type %in% c("clpm","dpm")){
    vals <- sapply(2:T, function(t){
      idx <- which(pe$lhs == paste0("y",t) & pe$rhs == paste0("x",t-1))
      if(length(idx) == 0) NA_real_ else pe$est[idx[1]]
    })
  } else { # riclpm: wy_t ~ wx_(t-1)
    vals <- sapply(2:T, function(t){
      idx <- which(pe$lhs == paste0("wy",t) & pe$rhs == paste0("wx",t-1))
      if(length(idx) == 0) NA_real_ else pe$est[idx[1]]
    })
  }

  if(all(is.na(vals))) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

############################################################
##  5. Data simulation with random nonlinearity + RW B_t
############################################################

simulate_panel_once <- function(
  N,
  T,
  k_lin,
  ax, ay, bx, by, rho,
  s0, eta0,
  rw_sd,
  shrink,
  nonlin_type = c("cubic","exp","sin","tanh"),
  exp_rate  = 0.7,
  sin_omega = 3,
  tanh_a    = 2
){
  nonlin_type <- match.arg(nonlin_type)

  # autoregressive matrix A
  A <- matrix(c(ax, by, bx, ay), 2, 2, byrow = TRUE)

  # 1) Base linear confounders
  Psi_base <- diag(rep(1, k_lin))
  U_lin <- rmvnorm(N, rep(0,k_lin), Psi_base)
  colnames(U_lin) <- paste0("c",1:k_lin)

  # 2) Nonlinear components per type (used only inside simulation)
  build_nonlin <- function(col) {
    if(nonlin_type == "cubic"){
      q  <- scale(col^2, center=TRUE, scale=TRUE)
      c3 <- scale(col^3, center=TRUE, scale=TRUE)
      cbind(q, c3)
    } else if(nonlin_type == "exp"){
      z <- exp(exp_rate * col)
      z_perp <- resid(lm(z ~ col))
      scale(z_perp, center=TRUE, scale=TRUE)
    } else if(nonlin_type == "sin"){
      z <- sin(sin_omega * col)
      z_perp <- resid(lm(z ~ col))
      scale(z_perp, center=TRUE, scale=TRUE)
    } else if(nonlin_type == "tanh"){
      z <- tanh(tanh_a * col)
      z_perp <- resid(lm(z ~ col))
      scale(z_perp, center=TRUE, scale=TRUE)
    }
  }

  if(nonlin_type == "cubic"){
    U_nonlin <- do.call(cbind, lapply(1:k_lin, function(j){
      out <- build_nonlin(U_lin[,j])
      colnames(out) <- c(paste0("c",j,"_quad"), paste0("c",j,"_cubic"))
      out
    }))
    k_nonlin <- 2*k_lin
  } else {
    U_nonlin <- do.call(cbind, lapply(1:k_lin, function(j){
      out <- build_nonlin(U_lin[,j])
      colnames(out) <- paste0("c",j,"_",nonlin_type)
      out
    }))
    k_nonlin <- k_lin
  }

  U <- cbind(U_lin, U_nonlin)
  k_total <- ncol(U)
  Psi <- cov(U)

  # 3) Sample baseline B0 using s0 and eta0
  B0 <- sample_B_matrix(
    k_lin    = k_lin,
    k_nonlin = k_nonlin,
    s0       = s0,
    eta0     = eta0
  ) # 2 x (k_lin + k_nonlin)

  # safety: match dims to U
  if(ncol(B0) != k_total){
    if(ncol(B0) < k_total){
      B0 <- B0[, rep(1:ncol(B0), length.out = k_total), drop=FALSE]
    } else {
      B0 <- B0[, 1:k_total, drop=FALSE]
    }
  }

  # 4) Random walk Beta_list with shrinkage toward B0
  Beta_list <- vector("list", T)
  Beta_list[[1]] <- B0

  for(t in 2:T){
    step <- matrix(rnorm(2*k_total, mean=0, sd=rw_sd), 2, k_total)
    B_prev <- Beta_list[[t-1]]
    B_new  <- B_prev + step + shrink * (B0 - B_prev)
    Beta_list[[t]] <- B_new
  }

  # 5) Build S_U and solve for stationary covariance
  S_U_list <- lapply(Beta_list, function(Bt) Bt %*% Psi %*% t(Bt))
  S_U <- Reduce("+", S_U_list) / T

  find_c <- function(A, S_U, rho){
    f <- function(c){
      S_tar <- matrix(c(1,c,c,1),2,2)
      S_dyn <- S_tar - S_U
      Sig_e <- S_dyn - t(A) %*% S_dyn %*% A
      if(any(is.na(Sig_e))) return(NA_real_)
      corr <- Sig_e[1,2] / sqrt(Sig_e[1,1]*Sig_e[2,2])
      corr - rho
    }
    out <- tryCatch(
      uniroot(f, interval=c(-0.99,0.99))$root,
      error = function(e) NA_real_
    )
    out
  }

  c_stat <- find_c(A, S_U, rho)
  if(is.na(c_stat)) return(NULL)

  S_target <- matrix(c(1,c_stat,c_stat,1),2,2)
  S_dyn <- S_target - S_U
  S_dyn <- (S_dyn + t(S_dyn))/2

  # positive-definiteness check
  ev <- eigen(S_dyn, symmetric=TRUE)$values
  if(any(ev <= 1e-10)) return(NULL)

  Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
  Sigma_e <- (Sigma_e + t(Sigma_e))/2
  if(any(diag(Sigma_e) <= 0) || any(is.na(Sigma_e))) return(NULL)

  Sigma_e1 <- S_dyn
  Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

  # 6) Simulate panel data
  Ddyn <- tryCatch(
    rmvnorm(N, c(0,0), Sigma_e1),
    error = function(e) NULL
  )
  if(is.null(Ddyn)) return(NULL)

  df <- matrix(NA, nrow=N, ncol=2*T + k_total)
  colnames(df) <- c(
    paste0("x",1:T),
    paste0("y",1:T),
    colnames(U)
  )
  df[, (2*T+1):(2*T+k_total)] <- U

  # first wave
  obs1 <- Ddyn + U %*% t(Beta_list[[1]])
  df[,1]   <- obs1[,1]
  df[,1+T] <- obs1[,2]

  for(t in 2:T){
    eps_t <- tryCatch(
      rmvnorm(N, c(0,0), Sigma_e),
      error = function(e) NULL
    )
    if(is.null(eps_t)) return(NULL)

    Ddyn <- Ddyn %*% t(A) + eps_t
    obs  <- Ddyn + U %*% t(Beta_list[[t]])
    df[,t]     <- obs[,1]
    df[,t + T] <- obs[,2]
  }

  list(
    df = as.data.frame(df),
    B_list = Beta_list,
    nl_type = nonlin_type
  )
}

############################################################
##  6. One replication: simulate + fit 6 models + extract CLs
############################################################

run_one_rep <- function(
  rep_id,
  N, T,
  k_lin,
  ax, ay, bx, by, rho,
  s0, eta0,
  rw_sd, shrink,
  model_clpm,
  model_riclpm,
  model_dpm,
  base_seed = 12345
){
  set.seed(base_seed + rep_id)

  # randomly pick nonlinearity (never purely linear)
  nonlin_type <- sample(c("cubic","exp","sin","tanh"), size=1)

  sim <- simulate_panel_once(
    N = N,
    T = T,
    k_lin = k_lin,
    ax = ax, ay = ay, bx = bx, by = by, rho = rho,
    s0 = s0,
    eta0 = eta0,
    rw_sd = rw_sd, shrink = shrink,
    nonlin_type = nonlin_type
  )

  if(is.null(sim)){
    return(data.frame(
      rep            = rep_id,
      nl_type        = nonlin_type,
      clpm_raw_est   = NA_real_,
      riclpm_est     = NA_real_,
      dpm_est        = NA_real_,
      clpm_lin_est   = NA_real_,
      clpm_xgb_est   = NA_real_,
      clpm_lasso_est = NA_real_,
      sd_B           = NA_real_,
      mean_B         = NA_real_,
      is_na_run      = 1
    ))
  }

  df      <- sim$df
  B_list  <- sim$B_list

  # compute sd_B and mean_B across time and entries
  sd_B <- NA_real_
  mean_B <- NA_real_
  if(!is.null(B_list) && length(B_list) == T){
    B_array <- tryCatch(
      simplify2array(B_list),
      error=function(e) NULL
    )
    if(!is.null(B_array)){
      # B_array: 2 x k x T
      sd_array <- apply(B_array, c(1,2), sd, na.rm=TRUE)
      sd_B <- mean(sd_array, na.rm=TRUE)
      mean_B <- mean(B_array, na.rm=TRUE)
    }
  }

  # CLPM raw
  fit_clpm_raw <- tryCatch(
    fit_clpm(df, model_clpm),
    error = function(e) NULL
  )
  clpm_raw_est <- if(is.null(fit_clpm_raw)) NA_real_ else
    tryCatch(extract_crosslag_mean(fit_clpm_raw, T, "clpm"),
             error = function(e) NA_real_)

  # RI-CLPM
  fit_ric <- tryCatch(
    fit_riclpm(df, model_riclpm),
    error = function(e) NULL
  )
  riclpm_est <- if(is.null(fit_ric)) NA_real_ else
    tryCatch(extract_crosslag_mean(fit_ric, T, "riclpm"),
             error = function(e) NA_real_)

  # DPM
  fit_d <- tryCatch(
    fit_dpm(df, model_dpm),
    error = function(e) NULL
  )
  dpm_est <- if(is.null(fit_d)) NA_real_ else
    tryCatch(extract_crosslag_mean(fit_d, T, "dpm"),
             error = function(e) NA_real_)

  # CLPM after linear residualisation
  df_lin <- tryCatch(
    residualise_panel_linearC(df),
    error = function(e) NULL
  )
  clpm_lin_est <- if(is.null(df_lin)) NA_real_ else {
    fit <- tryCatch(fit_clpm(df_lin, model_clpm), error=function(e) NULL)
    if(is.null(fit)) NA_real_ else
      tryCatch(extract_crosslag_mean(fit, T, "clpm"),
               error = function(e) NA_real_)
  }

  # CLPM after XGBoost residualisation
  df_xgb <- tryCatch(
    residualise_panel_xgb(df),
    error = function(e) NULL
  )
  clpm_xgb_est <- if(is.null(df_xgb)) NA_real_ else {
    fit <- tryCatch(fit_clpm(df_xgb, model_clpm), error=function(e) NULL)
    if(is.null(fit)) NA_real_ else
      tryCatch(extract_crosslag_mean(fit, T, "clpm"),
               error = function(e) NA_real_)
  }

  # CLPM after LASSO residualisation
  df_lasso <- tryCatch(
    residualise_panel_lasso(df),
    error = function(e) NULL
  )
  clpm_lasso_est <- if(is.null(df_lasso)) NA_real_ else {
    fit <- tryCatch(fit_clpm(df_lasso, model_clpm), error=function(e) NULL)
    if(is.null(fit)) NA_real_ else
      tryCatch(extract_crosslag_mean(fit, T, "clpm"),
               error = function(e) NA_real_)
  }

  data.frame(
    rep            = rep_id,
    nl_type        = nonlin_type,
    clpm_raw_est   = clpm_raw_est,
    riclpm_est     = riclpm_est,
    dpm_est        = dpm_est,
    clpm_lin_est   = clpm_lin_est,
    clpm_xgb_est   = clpm_xgb_est,
    clpm_lasso_est = clpm_lasso_est,
    sd_B           = sd_B,
    mean_B         = mean_B,
    is_na_run      = as.integer(all(is.na(c(
      clpm_raw_est, riclpm_est, dpm_est,
      clpm_lin_est, clpm_xgb_est, clpm_lasso_est
    ))))
  )
}

############################################################
##  7. Main wrapper: run many reps in parallel and summarise
############################################################

run_simulation_bias <- function(
  reps      = 200,
  N         = 10000,
  T         = 5,
  k_lin     = 3,
  ax        = 0.25,
  ay        = 0.25,
  bx        = 0.10,   # true X<-Y or Y<-X cross-lag
  by        = 0.10,   # true X->Y cross-lag (used as "true" in bias)
  rho       = 0.30,
  s0        = 0.20,   # total confounder strength
  eta0      = 0.30,   # nonlinearity fraction
  rw_sd     = 0.03,   # RW step SD for B_t
  shrink    = 0.10,   # shrinkage toward B0
  cores     = max(1, detectCores() - 1),
  base_seed = 12345
){
  start_time <- Sys.time()

  # model strings (once)
  model_clpm   <- build_clpm_model(T)
  model_riclpm <- build_riclpm_model(T)
  model_dpm    <- build_dpm_model(T)

  # parallel setup
  tasks <- 1:reps
  total_tasks <- length(tasks)

  cl <- makeCluster(cores)
  clusterExport(cl, c(
    # functions
    "sample_B_matrix",
    "residualise_panel_linearC",
    "residualise_panel_xgb",
    "residualise_panel_lasso",
    "build_clpm_model",
    "build_riclpm_model",
    "build_dpm_model",
    "fit_clpm",
    "fit_riclpm",
    "fit_dpm",
    "extract_crosslag_mean",
    "simulate_panel_once",
    "run_one_rep",
    # parameters
    "N","T","k_lin",
    "ax","ay","bx","by","rho",
    "s0","eta0","rw_sd","shrink",
    "model_clpm","model_riclpm","model_dpm",
    "base_seed"
  ), envir=environment())

  clusterEvalQ(cl, {
    library(mvtnorm)
    library(lavaan)
    library(xgboost)
    library(glmnet)
  })

  pb <- txtProgressBar(min=0, max=total_tasks, style=3)
  progress_counter <- 0

  idx <- 1:total_tasks
  chunk_size <- max(1, floor(total_tasks / 50))
  chunks <- split(idx, ceiling(seq_along(idx) / chunk_size))

  results_list <- vector("list", length(chunks))

  for(ci in seq_along(chunks)){
    idx_vec <- chunks[[ci]]

    res_chunk <- parLapplyLB(
      cl,
      idx_vec,
      function(i){
        run_one_rep(
          rep_id = i,
          N = N, T = T, k_lin = k_lin,
          ax = ax, ay = ay, bx = bx, by = by, rho = rho,
          s0 = s0, eta0 = eta0,
          rw_sd = rw_sd, shrink = shrink,
          model_clpm = model_clpm,
          model_riclpm = model_riclpm,
          model_dpm = model_dpm,
          base_seed = base_seed
        )
      }
    )

    results_list[[ci]] <- res_chunk

    progress_counter <- progress_counter + length(idx_vec)
    setTxtProgressBar(pb, progress_counter)
  }

  close(pb)
  stopCluster(cl)

  raw_df <- bind_rows(unlist(results_list, recursive = FALSE))

  # summary: mean estimate and bias per model
  summary_df <- raw_df %>%
    mutate(
      bias_clpm_raw   = clpm_raw_est   - by,
      bias_riclpm     = riclpm_est     - by,
      bias_dpm        = dpm_est        - by,
      bias_clpm_lin   = clpm_lin_est   - by,
      bias_clpm_xgb   = clpm_xgb_est   - by,
      bias_clpm_lasso = clpm_lasso_est - by
    ) %>%
    summarise(
      mean_est_clpm_raw   = mean(clpm_raw_est,   na.rm=TRUE),
      mean_est_riclpm     = mean(riclpm_est,     na.rm=TRUE),
      mean_est_dpm        = mean(dpm_est,        na.rm=TRUE),
      mean_est_clpm_lin   = mean(clpm_lin_est,   na.rm=TRUE),
      mean_est_clpm_xgb   = mean(clpm_xgb_est,   na.rm=TRUE),
      mean_est_clpm_lasso = mean(clpm_lasso_est, na.rm=TRUE),

      mean_bias_clpm_raw   = mean(bias_clpm_raw,   na.rm=TRUE),
      mean_bias_riclpm     = mean(bias_riclpm,     na.rm=TRUE),
      mean_bias_dpm        = mean(bias_dpm,        na.rm=TRUE),
      mean_bias_clpm_lin   = mean(bias_clpm_lin,   na.rm=TRUE),
      mean_bias_clpm_xgb   = mean(bias_clpm_xgb,   na.rm=TRUE),
      mean_bias_clpm_lasso = mean(bias_clpm_lasso, na.rm=TRUE),

      mean_absbias_clpm_raw   = mean(abs(bias_clpm_raw),   na.rm=TRUE),
      mean_absbias_riclpm     = mean(abs(bias_riclpm),     na.rm=TRUE),
      mean_absbias_dpm        = mean(abs(bias_dpm),        na.rm=TRUE),
      mean_absbias_clpm_lin   = mean(abs(bias_clpm_lin),   na.rm=TRUE),
      mean_absbias_clpm_xgb   = mean(abs(bias_clpm_xgb),   na.rm=TRUE),
      mean_absbias_clpm_lasso = mean(abs(bias_clpm_lasso), na.rm=TRUE),

      sd_est_clpm_raw   = sd(clpm_raw_est,   na.rm=TRUE),
      sd_est_riclpm     = sd(riclpm_est,     na.rm=TRUE),
      sd_est_dpm        = sd(dpm_est,        na.rm=TRUE),
      sd_est_clpm_lin   = sd(clpm_lin_est,   na.rm=TRUE),
      sd_est_clpm_xgb   = sd(clpm_xgb_est,   na.rm=TRUE),
      sd_est_clpm_lasso = sd(clpm_lasso_est, na.rm=TRUE),

      mean_sd_B   = mean(sd_B,   na.rm=TRUE),
      mean_mean_B = mean(mean_B, na.rm=TRUE),

      na_runs    = sum(is_na_run),
      total_runs = n(),
      na_rate    = na_runs / total_runs
    )

  runtime <- Sys.time() - start_time

  list(
    raw      = raw_df,
    summary  = summary_df,
    runtime  = runtime
  )
}

############################################################
##  8. Plotting helpers (general + per-type)
############################################################

# Helper: general long-format conversion for summary
.make_plot_df <- function(summary, true_by){
  data.frame(
    method = factor(c(
      "CLPM (raw)",
      "RI-CLPM",
      "DPM",
      "CLPM + Linear Residualisation",
      "CLPM + XGBoost Residualisation",
      "CLPM + LASSO Residualisation"
    ), levels = c(
      "CLPM (raw)",
      "RI-CLPM",
      "DPM",
      "CLPM + Linear Residualisation",
      "CLPM + XGBoost Residualisation",
      "CLPM + LASSO Residualisation"
    )),
    mean_bias = c(
      summary$mean_bias_clpm_raw,
      summary$mean_bias_riclpm,
      summary$mean_bias_dpm,
      summary$mean_bias_clpm_lin,
      summary$mean_bias_clpm_xgb,
      summary$mean_bias_clpm_lasso
    ),
    sd_est = c(
      summary$sd_est_clpm_raw,
      summary$sd_est_riclpm,
      summary$sd_est_dpm,
      summary$sd_est_clpm_lin,
      summary$sd_est_clpm_xgb,
      summary$sd_est_clpm_lasso
    ),
    mean_abs_bias = c(
      summary$mean_absbias_clpm_raw,
      summary$mean_absbias_riclpm,
      summary$mean_absbias_dpm,
      summary$mean_absbias_clpm_lin,
      summary$mean_absbias_clpm_xgb,
      summary$mean_absbias_clpm_lasso
    ),
    mean_rel_bias = c(
      summary$mean_bias_clpm_raw   / true_by,
      summary$mean_bias_riclpm     / true_by,
      summary$mean_bias_dpm        / true_by,
      summary$mean_bias_clpm_lin   / true_by,
      summary$mean_bias_clpm_xgb   / true_by,
      summary$mean_bias_clpm_lasso / true_by
    )
  )
}

# 8.1 General: mean bias + SD
plot_bias_results <- function(sim_out, true_by){

  df_plot <- .make_plot_df(sim_out$summary, true_by)

  ggplot(df_plot, aes(x = method, y = mean_bias)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean_bias - sd_est,
                      ymax = mean_bias + sd_est),
                  width = 0.15, size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title="Mean Bias with SD",
         subtitle=paste0("True β = ", true_by),
         y="Bias", x="Method") +
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

# 8.2 General: mean absolute bias
plot_absbias_results <- function(sim_out, true_by){

  df_plot <- .make_plot_df(sim_out$summary, true_by)

  ggplot(df_plot, aes(x = method, y = mean_abs_bias)) +
    geom_point(size = 4) +
    labs(title="Mean Absolute Bias",
         subtitle=paste0("True β = ", true_by),
         y="Absolute Bias", x="Method") +
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

# 8.3 General: mean relative bias
plot_relbias_results <- function(sim_out, true_by){

  df_plot <- .make_plot_df(sim_out$summary, true_by)

  ggplot(df_plot, aes(x = method, y = mean_rel_bias)) +
    geom_point(size = 4, color="purple") +
    geom_hline(yintercept = 0, color = "red", linetype="dashed") +
    labs(title="Mean Relative Bias",
         subtitle=paste0("True β = ", true_by),
         y="Relative Bias (bias / true β)", x="Method") +
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

############################################################
##  9. Per–nonlinearity-type summary and plots
############################################################

make_type_summary <- function(raw_df, true_by){

  raw_df %>%
    group_by(nl_type) %>%
    summarise(
      mean_bias_clpm_raw   = mean(clpm_raw_est   - true_by, na.rm=TRUE),
      sd_est_clpm_raw      = sd(clpm_raw_est,           na.rm=TRUE),
      absbias_clpm_raw     = mean(abs(clpm_raw_est   - true_by), na.rm=TRUE),
      relbias_clpm_raw     = mean((clpm_raw_est   - true_by)/true_by, na.rm=TRUE),

      mean_bias_riclpm     = mean(riclpm_est     - true_by, na.rm=TRUE),
      sd_est_riclpm        = sd(riclpm_est,             na.rm=TRUE),
      absbias_riclpm       = mean(abs(riclpm_est     - true_by), na.rm=TRUE),
      relbias_riclpm       = mean((riclpm_est     - true_by)/true_by, na.rm=TRUE),

      mean_bias_dpm        = mean(dpm_est        - true_by, na.rm=TRUE),
      sd_est_dpm           = sd(dpm_est,                na.rm=TRUE),
      absbias_dpm          = mean(abs(dpm_est        - true_by), na.rm=TRUE),
      relbias_dpm          = mean((dpm_est        - true_by)/true_by, na.rm=TRUE),

      mean_bias_clpm_lin   = mean(clpm_lin_est   - true_by, na.rm=TRUE),
      sd_est_clpm_lin      = sd(clpm_lin_est,           na.rm=TRUE),
      absbias_clpm_lin     = mean(abs(clpm_lin_est   - true_by), na.rm=TRUE),
      relbias_clpm_lin     = mean((clpm_lin_est   - true_by)/true_by, na.rm=TRUE),

      mean_bias_clpm_xgb   = mean(clpm_xgb_est   - true_by, na.rm=TRUE),
      sd_est_clpm_xgb      = sd(clpm_xgb_est,           na.rm=TRUE),
      absbias_clpm_xgb     = mean(abs(clpm_xgb_est   - true_by), na.rm=TRUE),
      relbias_clpm_xgb     = mean((clpm_xgb_est   - true_by)/true_by, na.rm=TRUE),

      mean_bias_clpm_lasso = mean(clpm_lasso_est - true_by, na.rm=TRUE),
      sd_est_clpm_lasso    = sd(clpm_lasso_est,         na.rm=TRUE),
      absbias_clpm_lasso   = mean(abs(clpm_lasso_est - true_by), na.rm=TRUE),
      relbias_clpm_lasso   = mean((clpm_lasso_est - true_by)/true_by, na.rm=TRUE)
    )
}

.make_type_plot_df <- function(S_row){
  data.frame(
    nl_type = S_row$nl_type,
    method = factor(c(
      "CLPM (raw)","RI-CLPM","DPM",
      "CLPM + Linear","CLPM + XGBoost","CLPM + LASSO"
    ), levels = c(
      "CLPM (raw)","RI-CLPM","DPM",
      "CLPM + Linear","CLPM + XGBoost","CLPM + LASSO"
    )),
    mean_bias = c(
      S_row$mean_bias_clpm_raw,
      S_row$mean_bias_riclpm,
      S_row$mean_bias_dpm,
      S_row$mean_bias_clpm_lin,
      S_row$mean_bias_clpm_xgb,
      S_row$mean_bias_clpm_lasso
    ),
    sd_est = c(
      S_row$sd_est_clpm_raw,
      S_row$sd_est_riclpm,
      S_row$sd_est_dpm,
      S_row$sd_est_clpm_lin,
      S_row$sd_est_clpm_xgb,
      S_row$sd_est_clpm_lasso
    ),
    abs_bias = c(
      S_row$absbias_clpm_raw,
      S_row$absbias_riclpm,
      S_row$absbias_dpm,
      S_row$absbias_clpm_lin,
      S_row$absbias_clpm_xgb,
      S_row$absbias_clpm_lasso
    ),
    rel_bias = c(
      S_row$relbias_clpm_raw,
      S_row$relbias_riclpm,
      S_row$relbias_dpm,
      S_row$relbias_clpm_lin,
      S_row$relbias_clpm_xgb,
      S_row$relbias_clpm_lasso
    )
  )
}

plot_bias_by_type <- function(sim_out, true_by){

  S <- make_type_summary(sim_out$raw, true_by)

  df <- bind_rows(lapply(1:nrow(S), function(i) .make_type_plot_df(S[i,])))

  ggplot(df, aes(x=method, y=mean_bias, color=nl_type)) +
    geom_point(size=3) +
    geom_errorbar(aes(ymin = mean_bias - sd_est,
                      ymax = mean_bias + sd_est),
                  width = 0.15) +
    facet_wrap(~nl_type, ncol=2) +
    geom_hline(yintercept = 0, color="red", linetype="dashed") +
    labs(title="Bias by Nonlinearity Type",
         subtitle=paste0("True β = ", true_by),
         y="Bias", x="Method") +
    theme_minimal(base_size=13) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

plot_absbias_by_type <- function(sim_out, true_by){

  S <- make_type_summary(sim_out$raw, true_by)
  df <- bind_rows(lapply(1:nrow(S), function(i) .make_type_plot_df(S[i,])))

  ggplot(df, aes(x=method, y=abs_bias, color=nl_type)) +
    geom_point(size=3) +
    facet_wrap(~nl_type, ncol=2) +
    labs(title="Absolute Bias by Nonlinearity Type",
         subtitle=paste0("True β = ", true_by),
         y="Absolute Bias", x="Method") +
    theme_minimal(base_size=13) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

plot_relbias_by_type <- function(sim_out, true_by){

  S <- make_type_summary(sim_out$raw, true_by)
  df <- bind_rows(lapply(1:nrow(S), function(i) .make_type_plot_df(S[i,])))

  ggplot(df, aes(x=method, y=rel_bias, color=nl_type)) +
    geom_point(size=3) +
    geom_hline(yintercept = 0, color="red", linetype="dashed") +
    facet_wrap(~nl_type, ncol=2) +
    labs(title="Relative Bias by Nonlinearity Type",
         subtitle=paste0("True β = ", true_by),
         y="Relative Bias (bias / true β)", x="Method") +
    theme_minimal(base_size=13) +
    theme(axis.text.x = element_text(angle=25, hjust=1))
}

############################################################
##  Example usage (commented)
############################################################
out <- run_simulation_bias(
   reps  = 1000,
   N     = 10000,
   T     = 5,
   k_lin = 3,
   ax    = 0.2,
   ay    = 0.2,
   bx    = 0.10,
   by    = 0.10,   
   rho   = 0.1,
   s0    = 0.25,
   eta0  = 0.40,
   rw_sd = 0.07,
   shrink = 0.15,
   cores  = 7
 )
out$summary
out$summary$na_runs   
out$summary$runtime
plot_bias_results(out, true_by = 0.10)
plot_absbias_results(out, true_by = 0.10)
plot_relbias_results(out, true_by = 0.10)
plot_bias_by_type(out, true_by = 0.10)
plot_absbias_by_type(out, true_by = 0.10)
plot_relbias_by_type(out, true_by = 0.10)
