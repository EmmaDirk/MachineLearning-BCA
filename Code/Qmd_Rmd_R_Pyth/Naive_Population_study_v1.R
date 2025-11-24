run_one_condition <- function(
  N    = 2000,
  T    = 5,
  s0   = 0.15,   # total confounder strength
  eta0 = 0.30,   # non-linearity index
  k_lin = 3,
  ax   = 0.25,
  ay   = 0.25,
  bx   = 0.10,
  by   = 0.10,
  rho  = 0.30
) {

  ## ---- libs ----
  library(mvtnorm)
  library(lavaan)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(parallel)

  ## ---- CLPM dynamic matrix ----
  A <- matrix(c(ax, by,
                bx, ay), nrow = 2, byrow = TRUE)

  k_nonlin <- k_lin  # for now: one nonlinear component per confounder
  k_total  <- k_lin + k_nonlin

  ## ---------- helper: sample gamma_x, gamma_y ----------
  sample_gamma <- function(
      k_lin, k_nonlin, s0, eta0,
      min_abs = 0.01, max_abs = 0.30,
      max_tries = 1e6
  ){
    LX <- (1 - eta0) * s0 / 2
    LY <- (1 - eta0) * s0 / 2
    NX <- eta0 * s0 / 2
    NY <- eta0 * s0 / 2

    for (try in 1:max_tries) {
      u_x_lin  <- rnorm(k_lin)
      u_y_lin  <- rnorm(k_lin)
      u_x_nlin <- rnorm(k_nonlin)
      u_y_nlin <- rnorm(k_nonlin)

      u_x_lin  <- u_x_lin  / sqrt(sum(u_x_lin^2))
      u_y_lin  <- u_y_lin  / sqrt(sum(u_y_lin^2))
      u_x_nlin <- u_x_nlin / sqrt(sum(u_x_nlin^2))
      u_y_nlin <- u_y_nlin / sqrt(sum(u_y_nlin^2))

      gamma_x <- c(sqrt(LX) * u_x_lin,
                   sqrt(NX) * u_x_nlin)
      gamma_y <- c(sqrt(LY) * u_y_lin,
                   sqrt(NY) * u_y_nlin)

      if (all(abs(gamma_x) >= min_abs & abs(gamma_x) <= max_abs &
              abs(gamma_y) >= min_abs & abs(gamma_y) <= max_abs)) {
        return(list(gamma_x = gamma_x, gamma_y = gamma_y))
      }
    }
    stop("Could not sample gamma within max_tries.")
  }

  ## ---------- helper: generic dynamic simulator given U and B ----------
  simulate_from_U_B <- function(U, B) {
    N  <- nrow(U)
    k  <- ncol(U)

    Psi <- cov(U)
    S_U <- B %*% Psi %*% t(B)   # 2x2 confounder variance

    find_c <- function(A, S_U, rho) {
      f <- function(c) {
        S_target_c <- matrix(c(1, c,
                               c, 1),
                             nrow = 2, byrow = TRUE)
        S_dyn_c <- S_target_c - S_U
        Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
        v1 <- Sigma_e_c[1,1]
        v2 <- Sigma_e_c[2,2]
        cov12 <- Sigma_e_c[1,2]
        corr_e <- cov12 / sqrt(v1 * v2)
        corr_e - rho
      }
      uniroot(f, interval = c(-0.99, 0.99))$root
    }

    c_stat <- find_c(A, S_U, rho)
    S_target <- matrix(c(1, c_stat,
                         c_stat, 1),
                       nrow = 2, byrow = TRUE)
    S_dyn <- S_target - S_U

    Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
    Sigma_e <- (Sigma_e + t(Sigma_e)) / 2

    Sigma_e1 <- (S_dyn + t(S_dyn)) / 2

    Ddyn <- rmvnorm(N, mean = c(0,0), sigma = Sigma_e1)

    df <- matrix(NA, nrow = N, ncol = 2*T + k)
    colnames(df) <- c(paste0("x", 1:T),
                      paste0("y", 1:T),
                      paste0("c", 1:k))
    df[, (2*T + 1):(2*T + k)] <- U

    obs1 <- Ddyn + U %*% t(B)
    df[,1]    <- obs1[,1]
    df[,1+T ] <- obs1[,2]

    for (ti in 2:T) {
      Ddyn <- Ddyn %*% t(A) + rmvnorm(N, mean = c(0,0), sigma = Sigma_e)
      obs  <- Ddyn + U %*% t(B)
      df[,ti]    <- obs[,1]
      df[,ti+T ] <- obs[,2]
    }

    as.data.frame(df)
  }

  ## ---------- six generators, all using simulate_from_U_B ----------

  # 1) pure linear confounding (only first k_lin gammas)
  gen_linear <- function(gamma_x, gamma_y) {
    U_lin <- rmvnorm(N, mean = rep(0, k_lin), sigma = diag(k_lin))
    B_lin <- rbind(gamma_x[1:k_lin], gamma_y[1:k_lin])
    simulate_from_U_B(U_lin, B_lin)
  }

  # 2) cubic
  gen_cubic <- function(gamma_x, gamma_y) {
    U_lin <- rmvnorm(N, rep(0, k_lin), diag(k_lin))

    U_quad <- apply(U_lin, 2, function(col) {
      z_raw  <- col^2
      z_perp <- resid(lm(z_raw ~ col))
      as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
    })

    U_cubic <- apply(U_lin, 2, function(col) {
      z_raw  <- col^3
      z_perp <- resid(lm(z_raw ~ col))
      as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
    })

    U <- cbind(U_lin, U_quad, U_cubic)
    B <- rbind(gamma_x, gamma_y)
    simulate_from_U_B(U, B)
  }

  # 3) cutoff / plateau
  gen_cutoff <- function(gamma_x, gamma_y) {
    U_lin <- rmvnorm(N, rep(0, k_lin), diag(k_lin))
    thr   <- 1
    U_cut <- pmax(U_lin - thr, 0)
    U <- cbind(U_lin, U_cut)
    B <- rbind(gamma_x, gamma_y)
    simulate_from_U_B(U, B)
  }

  # 4) exponential
  gen_exp <- function(gamma_x, gamma_y) {
    U_lin <- rmvnorm(N, rep(0, k_lin), diag(k_lin))
    b_exp <- 0.7
    U_exp <- apply(U_lin, 2, function(col) {
      z_raw  <- exp(b_exp * col)
      z_perp <- resid(lm(z_raw ~ col))
      as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
    })
    U <- cbind(U_lin, U_exp)
    B <- rbind(gamma_x, gamma_y)
    simulate_from_U_B(U, B)
  }

  # 5) sinusoid
  gen_sin <- function(gamma_x, gamma_y) {
    U_lin <- rmvnorm(N, rep(0, k_lin), diag(k_lin))
    omega <- 3
    U_sin <- apply(U_lin, 2, function(col) {
      z_raw  <- sin(omega * col)
      z_perp <- resid(lm(z_raw ~ col))
      as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
    })
    U <- cbind(U_lin, U_sin)
    B <- rbind(gamma_x, gamma_y)
    simulate_from_U_B(U, B)
  }

  # 6) tanh
  gen_tanh <- function(gamma_x, gamma_y) {
    U_lin  <- rmvnorm(N, rep(0, k_lin), diag(k_lin))
    a_tanh <- 2
    U_tanh <- apply(U_lin, 2, function(col) {
      z_raw  <- tanh(a_tanh * col)
      z_perp <- resid(lm(z_raw ~ col))
      as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
    })
    U <- cbind(U_lin, U_tanh)
    B <- rbind(gamma_x, gamma_y)
    simulate_from_U_B(U, B)
  }

  ## ---------- CLPM / RI-CLPM / DPM syntax ----------

  model_clpm <- "
  x2 + y2 ~ x1 + y1
  x3 + y3 ~ x2 + y2
  x4 + y4 ~ x3 + y3
  x5 + y5 ~ x4 + y4

  x1 ~~ y1
  x2 ~~ y2
  x3 ~~ y3
  x4 ~~ y4
  x5 ~~ y5

  x1 ~~ x1; y1 ~~ y1
  x2 ~~ x2; y2 ~~ y2
  x3 ~~ x3; y3 ~~ y3
  x4 ~~ x4; y4 ~~ y4
  x5 ~~ x5; y5 ~~ y5

  x1 ~ 1; x2 ~ 1; x3 ~ 1; x4 ~ 1; x5 ~ 1
  y1 ~ 1; y2 ~ 1; y3 ~ 1; y4 ~ 1; y5 ~ 1
  "

  model_riclpm <- "
  rix =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
  riy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5
  rix ~~ rix
  riy ~~ riy
  rix ~~ riy

  x1 ~~ 0*x1; x2 ~~ 0*x2; x3 ~~ 0*x3; x4 ~~ 0*x4; x5 ~~ 0*x5
  y1 ~~ 0*y1; y2 ~~ 0*y2; y3 ~~ 0*y3; y4 ~~ 0*y4; y5 ~~ 0*y5

  wx1 =~ 1*x1; wx2 =~ 1*x2; wx3 =~ 1*x3; wx4 =~ 1*x4; wx5 =~ 1*x5
  wy1 =~ 1*y1; wy2 =~ 1*y2; wy3 =~ 1*y3; wy4 =~ 1*y4; wy5 =~ 1*y5

  rix ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
  rix ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5
  riy ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
  riy ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5

  wx1 ~~ wx1; wx2 ~~ wx2; wx3 ~~ wx3; wx4 ~~ wx4; wx5 ~~ wx5
  wy1 ~~ wy1; wy2 ~~ wy2; wy3 ~~ wy3; wy4 ~~ wy4; wy5 ~~ wy5
  wy1 ~~ wx1; wy2 ~~ wx2; wy3 ~~ wx3; wy4 ~~ wx4; wy5 ~~ wx5

  wx2 ~ wx1 + wy1
  wy2 ~ wx1 + wy1
  wx3 ~ wx2 + wy2
  wy3 ~ wx2 + wy2
  wx4 ~ wx3 + wy3
  wy4 ~ wx3 + wy3
  wx5 ~ wx4 + wy4
  wy5 ~ wx4 + wy4

  x1 + x2 + x3 + x4 + x5 ~ mx*1
  y1 + y2 + y3 + y4 + y5 ~ my*1
  "

  model_dpm <- "
  fx =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
  fy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5

  fx ~~ x1 + y1
  fy ~~ x1 + y1

  x2 + y2 ~ x1 + y1
  x3 + y3 ~ x2 + y2
  x4 + y4 ~ x3 + y3
  x5 + y5 ~ x4 + y4

  x1 ~~ y1
  x2 ~~ y2
  x3 ~~ y3
  x4 ~~ y4
  x5 ~~ y5

  fx ~~ fx
  fy ~~ fy
  fx ~~ fy

  x1 ~~ x1
  x2 ~~ x2
  x3 ~~ x3
  x4 ~~ x4
  x5 ~~ x5
  y1 ~~ y1
  y2 ~~ y2
  y3 ~~ y3
  y4 ~~ y4
  y5 ~~ y5

  x1 ~ 1
  x2 ~ 1
  x3 ~ 1
  x4 ~ 1
  x5 ~ 1
  y1 ~ 1
  y2 ~ 1
  y3 ~ 1
  y4 ~ 1
  y5 ~ 1
  "

  ## ---------- helper: BCA residualization ----------
  bca_residualize <- function(df, k_lin, method = c("lm", "lasso", "ranger", "xgb")) {
    method <- match.arg(method)
    conf_cols <- paste0("c", 1:k_lin)
    X_conf <- as.matrix(df[, conf_cols])

    get_resid <- function(y) {
      if (method == "lm") {
        fit <- lm(y ~ X_conf)
        return(resid(fit))
      } else if (method == "lasso") {
        cv <- glmnet::cv.glmnet(X_conf, y, alpha = 1, nfolds = 5)
        yhat <- predict(cv, newx = X_conf, s = "lambda.min")[,1]
        return(y - yhat)
      } else if (method == "ranger") {
        dat <- data.frame(y = y, X_conf)
        fit <- ranger::ranger(y ~ ., data = dat,
                              num.trees = 400,
                              mtry = max(1, floor(sqrt(ncol(X_conf)))),
                              min.node.size = 5)
        yhat <- predict(fit, dat)$predictions
        return(y - yhat)
      } else if (method == "xgb") {
        dtrain <- xgboost::xgb.DMatrix(data = X_conf, label = y)
        fit <- xgboost::xgb.train(
          params = list(
            objective = "reg:squarederror",
            eta       = 0.1,
            max_depth = 3
          ),
          data    = dtrain,
          nrounds = 200,
          verbose = 0
        )
        yhat <- predict(fit, X_conf)
        return(y - yhat)
      }
    }

    df2 <- df
    df2$x1 <- get_resid(df$x1)
    df2$y1 <- get_resid(df$y1)
    df2
  }

  ## ---------- helper: extract cross-lag estimates ----------
  get_xy_clpm <- function(fit) {
    pe <- parameterEstimates(fit)
    est <- pe$est[pe$lhs == "y2" & pe$rhs == "x1" & pe$op == "~"]
    ifelse(length(est) == 1, est, NA_real_)
  }

  get_xy_dpm <- function(fit) {
    pe <- parameterEstimates(fit)
    est <- pe$est[pe$lhs == "y2" & pe$rhs == "x1" & pe$op == "~"]
    ifelse(length(est) == 1, est, NA_real_)
  }

  get_xy_riclpm <- function(fit) {
    pe <- parameterEstimates(fit)
    est <- pe$est[pe$lhs == "wy2" & pe$rhs == "wx1" & pe$op == "~"]
    ifelse(length(est) == 1, est, NA_real_)
  }

  ## ---------- 1) sample gamma ----------
  g  <- sample_gamma(k_lin = k_lin, k_nonlin = k_nonlin, s0 = s0, eta0 = eta0)
  gx <- g$gamma_x
  gy <- g$gamma_y

  ## ---------- 2) simulate all 6 forms ----------
  df_lin  <- gen_linear(gx, gy)
  df_cub  <- gen_cubic(gx, gy)
  df_cut  <- gen_cutoff(gx, gy)
  df_exp  <- gen_exp(gx, gy)
  df_sin  <- gen_sin(gx, gy)
  df_tanh <- gen_tanh(gx, gy)

  sim_list <- list(
    linear = df_lin,
    cubic  = df_cub,
    cutoff = df_cut,
    exp    = df_exp,
    sin    = df_sin,
    tanh   = df_tanh
  )

  ## ---------- 3) fit all models per functional form ----------
  out_rows <- list()

  for (fn in names(sim_list)) {
    df <- sim_list[[fn]]

    # CLPM
    fit_clpm <- lavaan::sem(model_clpm, data = df, meanstructure = TRUE)
    est_clpm <- get_xy_clpm(fit_clpm)

    # RI-CLPM
    fit_riclpm <- lavaan::sem(model_riclpm, data = df, meanstructure = TRUE)
    est_riclpm <- get_xy_riclpm(fit_riclpm)

    # DPM
    fit_dpm <- lavaan::sem(model_dpm, data = df, meanstructure = TRUE)
    est_dpm <- get_xy_dpm(fit_dpm)

    # BCA variants
    df_lm    <- bca_residualize(df, k_lin, "lm")
    df_lasso <- bca_residualize(df, k_lin, "lasso")
    df_rf    <- bca_residualize(df, k_lin, "ranger")
    df_xgb   <- bca_residualize(df, k_lin, "xgb")

    fit_lm    <- lavaan::sem(model_clpm, data = df_lm,    meanstructure = TRUE)
    fit_lasso <- lavaan::sem(model_clpm, data = df_lasso, meanstructure = TRUE)
    fit_rf    <- lavaan::sem(model_clpm, data = df_rf,    meanstructure = TRUE)
    fit_xgb   <- lavaan::sem(model_clpm, data = df_xgb,   meanstructure = TRUE)

    out_rows[[fn]] <- data.frame(
      form          = fn,
      est_clpm      = est_clpm,
      est_riclpm    = est_riclpm,
      est_dpm       = est_dpm,
      est_bca_lm    = get_xy_clpm(fit_lm),
      est_bca_lasso = get_xy_clpm(fit_lasso),
      est_bca_ranger= get_xy_clpm(fit_rf),
      est_bca_xgb   = get_xy_clpm(fit_xgb),
      N      = N,
      T      = T,
      k_lin  = k_lin,
      s0     = s0,
      eta0   = eta0,
      ax     = ax,
      ay     = ay,
      bx     = bx,
      by     = by,
      rho    = rho,
      gamma_x = I(list(gx)),
      gamma_y = I(list(gy)),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_rows)
}

cl <- makeCluster(7)

clusterEvalQ(cl, {
  library(mvtnorm)
  library(lavaan)
  library(glmnet)
  library(ranger)
  library(xgboost)
})

clusterExport(cl, "run_one_condition")

res_list <- parLapply(cl, 1:100, function(r) {
  set.seed(1000 + r)
  run_one_condition(
    N    = 2000,
    T    = 5,
    s0   = 0.15,
    eta0 = 0.30,
    k_lin = 3
  )
})

stopCluster(cl)

res_all <- do.call(rbind, res_list)
