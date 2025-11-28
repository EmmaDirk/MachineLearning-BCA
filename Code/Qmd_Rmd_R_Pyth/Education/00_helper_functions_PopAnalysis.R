##############################################################################--
##### helper functions


# compute correlation matrix between continuous and binary variables -----------
# split at 0
transform_mixed_corr <- function(Sigma, n_cont = 50, n_bin = 50) {
  stopifnot(is.matrix(Sigma), nrow(Sigma) == ncol(Sigma))
  p <- nrow(Sigma)
  stopifnot(n_cont + n_bin == p)
  # indices
  I <- seq_len(n_cont)
  J <- (n_cont + 1):p

  # blocks of latent correlations
  CC <- Sigma[I, I, drop = FALSE]
  CB <- Sigma[I, J, drop = FALSE]
  BB <- Sigma[J, J, drop = FALSE]

  # constants
  c_cb <- sqrt(2 / pi)     # polyserial factor at t = 0
  # guard against tiny numeric excursions of |rho|>1
  clamp <- function(x) pmax(pmin(x, 1), -1)

  # transform blocks
  R_CC <- CC
  R_CB <- c_cb * CB
  R_BC <- t(R_CB)
  R_BB <- (2 / pi) * asin(clamp(BB))

  # assemble and fix diagonals
  R_obs <- matrix(NA_real_, p, p)
  R_obs[I, I] <- R_CC
  R_obs[I, J] <- R_CB
  R_obs[J, I] <- R_BC
  R_obs[J, J] <- R_BB

  diag(R_obs) <- 1
  # enforce symmetry (defensive against tiny asymmetries)
  R_obs <- (R_obs + t(R_obs)) / 2

  R_obs
}


# compute coefficients cross lagged structure ----------------------------------
computeSigmaCoefficients_linear <- function(
    nT,                # number of post-baseline waves
    Sigma_C,           # covariance of baseline covariates C (no expansion)
    BetaC0,            # pX x pC regression at time 0 (rows = X,Y; cols = C)
    BetaC_List,        # list of length nT with per-wave BetaC matrices (pX x pC)
    CrossLagged_List,  # list of length nT with per-wave cross-lag matrices (pX x pX)
    SigmaTarget        # target Var(X_k) under stationarity at each wave (pX x pX)
){
  # pre-processing (allow single-matrix inputs replicated across time)
  if(!is.list(BetaC_List))        BetaC_List        <- list(BetaC_List)
  if(!is.list(CrossLagged_List))  CrossLagged_List  <- list(CrossLagged_List)
  if(length(BetaC_List) == 1)        BetaC_List     <- replicate(nT, list(BetaC_List[[1]]))
  if(length(CrossLagged_List) == 1)  CrossLagged_List <- replicate(nT, list(CrossLagged_List[[1]]))

  # No nonlinear expansion:
  Sigma_C_lin <- Sigma_C
  pC <- ncol(Sigma_C_lin)
  colnames(Sigma_C_lin) <- rownames(Sigma_C_lin) <- paste0("C", 1:pC)
  VarExpC_List <- lapply(BetaC_List, FUN = function(B) B %*% Sigma_C_lin %*% t(B))

  Tm <- nT + 1    # include time 0
  # 1) Cov[X, C] for each wave k=0…nT ------------------------------------------
  Cov_XC <- vector("list", Tm)
  # time 0
  Cov_XC[[1]] <- BetaC0 %*% Sigma_C_lin

  # subsequent times 1…nT
  for(k in 2:Tm) {
    # Cov(X_k, C) = A_k-1 Cov(X_{k-1}, C) + B_k-1 Sigma_C
    Cov_XC[[k]] <-
      CrossLagged_List[[k-1]] %*% Cov_XC[[k-1]] +
      BetaC_List[[k-1]]       %*% Sigma_C_lin
  }

  # 2) Cov[X, X] between every pair of waves -----------------------------------
  Cov_XX <- vector("list", Tm)
  for(i in 1:Tm) Cov_XX[[i]] <- vector("list", Tm)

  # variance at time 0
  VarExpl0     <- BetaC0 %*% Sigma_C_lin %*% t(BetaC0)
  SigmaEps_T0  <- SigmaTarget - VarExpl0  # enforce Var(X_0)=SigmaTarget
  Cov_XX[[1]][[1]] <- VarExpl0 + SigmaEps_T0

  Sigma_Eps_List <- vector("list", nT)

  # times 1…nT
  for(k in 2:Tm) {
    A <- CrossLagged_List[[k-1]]
    B <- BetaC_List[[k-1]]

    # 2a) Var[X_k]
    VarExpl <- A %*% Cov_XX[[k-1]][[k-1]] %*% t(A) +   # propagated Var
      A %*% Cov_XC[[k-1]]         %*% t(B) +  # cross-term #1
      B %*% t(Cov_XC[[k-1]])      %*% t(A) +  # cross-term #2
      B %*% Sigma_C_lin           %*% t(B)    # new variance from C

    SigmaEps_T <- SigmaTarget - VarExpl
    Sigma_Eps_List[[k-1]] <- SigmaEps_T
    Cov_XX[[k]][[k]] <- VarExpl + SigmaEps_T

    # 2b) Cov[X_k, X_j] for j<k
    for(j in 1:(k-1)) {
      Cov_XX[[k]][[j]] <-
        A %*% Cov_XX[[k-1]][[j]] +   # propagate old covariances
        B %*% t(Cov_XC[[j]])         # contemporaneous regressor cross-term
      Cov_XX[[j]][[k]] <- t(Cov_XX[[k]][[j]])
    }
  }

  # 3) assemble full Cov((C, X0…XnT))
  pC <- ncol(Sigma_C_lin)
  pX <- nrow(CrossLagged_List[[1]])
  P  <- pC + pX*Tm
  Cov_full <- matrix(0, P, P)

  # CC block
  Cov_full[1:pC, 1:pC] <- Sigma_C_lin

  # CX and XC blocks
  for(k in 1:Tm) {
    rowsX <- pC + ((k-1)*pX + 1:pX)
    Cov_full[rowsX, 1:pC] <- Cov_XC[[k]]
    Cov_full[1:pC, rowsX] <- t(Cov_XC[[k]])
  }

  # XX blocks
  for(k in 1:Tm) for(j in 1:Tm) {
    rowsK <- pC + ((k-1)*pX + 1:pX)
    colsJ <- pC + ((j-1)*pX + 1:pX)
    Cov_full[rowsK, colsJ] <- Cov_XX[[k]][[j]]
  }

  colnames(Cov_full) <- rownames(Cov_full) <- c(
    colnames(Sigma_C_lin),
    paste0(c("X","Y"), rep(0:nT, each = 2))
  )

  list(
    Sigma_full        = Cov_full,
    CrossLagged_List  = CrossLagged_List,
    Sigma_Eps_List    = Sigma_Eps_List,
    BetaC_List        = BetaC_List,
    BetaC0            = BetaC0,
    Sigma_C           = Sigma_C_lin,
    SigmaTarget       = SigmaTarget,
    VarExpC_List      = VarExpC_List
  )
}

computeSigmaWithBurnin_linear <- function(
    nT,
    burnin = 50,
    BetaC_List,
    CrossLagged_List,
    Sigma_C,
    SigmaTarget
){
  # pre-processing (allow single-matrix inputs replicated across time)
  if(!is.list(BetaC_List))        BetaC_List        <- list(BetaC_List)
  if(!is.list(CrossLagged_List))  CrossLagged_List  <- list(CrossLagged_List)
  if(length(BetaC_List) == 1)        BetaC_List     <- replicate(nT, list(BetaC_List[[1]]))
  if(length(CrossLagged_List) == 1)  CrossLagged_List <- replicate(nT, list(CrossLagged_List[[1]]))

  T_real  <- length(BetaC_List) + 1
  nT_full <- burnin + T_real

  # Prepare extended series for burn-in + real waves
  Beta0      <- BetaC_List[[1]]
  p_baseline <- ncol(Beta0)
  CL0        <- CrossLagged_List[[1]]
  BetaC_ext  <- c(rep(list(Beta0), burnin+1), BetaC_List)
  CL_ext     <- c(rep(list(CL0),   burnin+1), CrossLagged_List)

  # Compute full covariance across all waves
  out_full <- computeSigmaCoefficients_linear(
    nT               = nT_full,
    Sigma_C          = Sigma_C,
    BetaC0           = Beta0,
    BetaC_List       = BetaC_ext,
    CrossLagged_List = CL_ext,
    SigmaTarget      = SigmaTarget
  )

  # Capture original baseline names
  if (is.null(colnames(out_full$Sigma_full))) {
    stop("Sigma_full must have column names for baseline extraction.")
  }
  baseline_names <- colnames(out_full$Sigma_full)[1:p_baseline]

  # Indices: keep baseline and post-burnin
  # Baseline occupy 1:p_baseline, waves each 2 vars
  start_post <- burnin * 2 + 1
  end_full   <- nT_full * 2
  idx_keep   <- c(seq_len(p_baseline), p_baseline + 2 + seq(start_post, end_full))

  # Trim the full covariance matrix
  out_full$Sigma_full <- out_full$Sigma_full[idx_keep, idx_keep]

  # Rename rows/cols: baseline + X0,Y0...X{T_real-1},Y{T_real-1}
  real_names <- paste0(rep(c("X","Y"), times = T_real),
                       rep(0:(T_real-1), each = 2))
  colnames(out_full$Sigma_full) <- rownames(out_full$Sigma_full) <-
    c(baseline_names, real_names)

  # Trim any per-wave lists
  if (!is.null(out_full$Sigma_Eps_List)) {
    out_full$Sigma_Eps_List <- out_full$Sigma_Eps_List[(burnin+1):nT_full]
  }
  if (!is.null(out_full$CrossLagged_List)) {
    out_full$CrossLagged_List <- out_full$CrossLagged_List[(burnin+1):nT_full]
  }
  if (!is.null(out_full$BetaC_List)) {
    out_full$BetaC_List <- out_full$BetaC_List[(burnin+1):nT_full]
  }
  if (!is.null(out_full$VarExpC_List)) {
    out_full$VarExpC_List <- out_full$VarExpC_List[(burnin+1):nT_full]
  }

  out_full
}




##### automate computation of Bs and compute ranks and cfa etc -----------------

automate_scenario_cfa <- function(
    BetaBase,                 # 2 x p matrix (rows = X,Y; cols = covariates)
    multi_List,     # List of multiplication factors
    corr_mat_theo,            # p x p correlation/covariance of covariates
    VarComp,                  # length 2*n_waves; diagonal variances to set
    p_dich = 50,            # optional: how many of the p covariates are dichotomous
    covar_names = NULL        # optional: names for covariates; length p
){
  # 2-factor CFA
  model2CFA <- "X =~ X1 + X2 + X3 + X4 + X5; Y =~ Y1 + Y2 + Y3 + Y4 + Y5"


  # ---- dependencies ----------------------------------------------------------
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.")
  if (!requireNamespace("lavaan", quietly = TRUE))
    stop("Package 'lavaan' is required.")

  # ---- sanity checks ---------------------------------------------------------
  if (!is.matrix(BetaBase) || nrow(BetaBase) != 2)
    stop("BetaBase must be a 2 x p matrix (rows correspond to X and Y).")
  p <- ncol(BetaBase)
  n_waves <- length(multi_List)

  if (!is.matrix(corr_mat_theo) || any(dim(corr_mat_theo) != c(p, p)))
    stop("corr_mat_theo must be a p x p matrix, with p = ncol(BetaBase).")

  # ---- build per-wave coefficient matrices -----------------------------------
  B_XY_List <- lapply(multi_List, FUN = function(M) M*BetaBase)


  # Stack (rbind) into a (2*n_waves) x p matrix
  B_mat <- do.call(rbind, B_XY_List)

  # ---- row and column names ---------------------------------------------------
  # Row names: X1,Y1,X2,Y2,...,Xn,Yn
  rn <- paste0(rep(c("X","Y"), times = n_waves),
               rep(seq_len(n_waves), each = 2))
  rownames(B_mat) <- rn

  # Column names: either provided or auto
  if (!is.null(covar_names)) {
    if (length(covar_names) != p) stop("covar_names must have length p.")
    colnames(B_mat) <- covar_names
  } else {
    # If p_dich is given, use your C / C_dich pattern; otherwise plain C1..Cp
    if (!is.null(p_dich)) {
      if (p_dich < 0 || p_dich > p) stop("p_dich must be between 0 and p.")
      base_names  <- paste0("C", seq_len(p))
      suffixes    <- c(rep("", p - p_dich), rep("_dich", p_dich))
      colnames(B_mat) <- paste0(base_names, suffixes)
    } else {
      colnames(B_mat) <- paste0("C", seq_len(p))
    }
  }

  # ---- sort rows by X first, then Y (your BX/BY step) ------------------------
  BX <- B_mat[grepl("^X", rownames(B_mat)), , drop = FALSE]
  BY <- B_mat[grepl("^Y", rownames(B_mat)), , drop = FALSE]
  B_sorted <- rbind(BX, BY)

  # ---- rank checks -----------------------------------------------------------
  rBX <- as.numeric(Matrix::rankMatrix(BX))
  rBY <- as.numeric(Matrix::rankMatrix(BY))
  rB  <- as.numeric(Matrix::rankMatrix(B_sorted))
  ranks <- c("Rank[BX]" = rBX, "Rank[BY]" = rBY, "Rank[BXY]" = rB)

  # ---- implied covariance and CFA fit ----------------------------------------
  # Expected covariance of (X1..Xn,Y1..Yn) explained by covariates:
  ExpVar <- B_sorted %*% corr_mat_theo %*% t(B_sorted)

  # Overwrite diagonal with VarComp (your sufficient condition construction)
  diag(ExpVar) <- VarComp

  # Fit the 2-factor CFA
  fit <- lavaan::sem(
    model = model2CFA,
    sample.cov = ExpVar,
    sample.cov.rescale = FALSE,
    sample.nobs = 1e5
  )

  fit_meas <- round(lavaan::fitMeasures(fit)[c("chisq", "pvalue",
                                               "cfi", "rmsea", "tli", "srmr")], 4)

  # ---- return ----------------------------------------------------------------
  out <- list(
    B_XY_List = B_XY_List,
    B_XY_Mat   = B_mat,
    B_sorted_X_then_Y = B_sorted,
    BX = BX, BY = BY,
    Ranks           = ranks,
    ExpVar          = ExpVar,
    lavaan_fit      = fit,
    fit_measures    = fit_meas
  )
  return(out)
}


# =============================================================================
# MODELS (observed variables CAPITALIZED)
# =============================================================================

model_clpm <- "
X2 + Y2 ~ X1 + Y1
X3 + Y3 ~ X2 + Y2
X4 + Y4 ~ X3 + Y3
X5 + Y5 ~ X4 + Y4

X1 ~~ Y1
X2 ~~ Y2
X3 ~~ Y3
X4 ~~ Y4
X5 ~~ Y5

X1 ~~ X1
Y1 ~~ Y1
X2 ~~ X2
Y2 ~~ Y2
X3 ~~ X3
Y3 ~~ Y3
X4 ~~ X4
Y4 ~~ Y4
X5 ~~ X5
Y5 ~~ Y5

X1 ~ 1
X2 ~ 1
X3 ~ 1
X4 ~ 1
X5 ~ 1
Y1 ~ 1
Y2 ~ 1
Y3 ~ 1
Y4 ~ 1
Y5 ~ 1
"

model_clpm_lag2 <- "
X2 + Y2 ~ X1 + Y1
X3 + Y3 ~ X2 + Y2 + X1 + Y1
X4 + Y4 ~ X3 + Y3 + X2 + Y2
X5 + Y5 ~ X4 + Y4 + X3 + Y3

X1 ~~ Y1
X2 ~~ Y2
X3 ~~ Y3
X4 ~~ Y4
X5 ~~ Y5

X1 ~~ X1
Y1 ~~ Y1
X2 ~~ X2
Y2 ~~ Y2
X3 ~~ X3
Y3 ~~ Y3
X4 ~~ X4
Y4 ~~ Y4
X5 ~~ X5
Y5 ~~ Y5

X1 ~ 1
X2 ~ 1
X3 ~ 1
X4 ~ 1
X5 ~ 1
Y1 ~ 1
Y2 ~ 1
Y3 ~ 1
Y4 ~ 1
Y5 ~ 1
"

model_dpm <- "
FX =~ 1*X2 + 1*X3 + 1*X4 + 1*X5
FY =~ 1*Y2 + 1*Y3 + 1*Y4 + 1*Y5

FX ~~ X1 + Y1
FY ~~ X1 + Y1

X2 + Y2 ~ X1 + Y1
X3 + Y3 ~ X2 + Y2
X4 + Y4 ~ X3 + Y3
X5 + Y5 ~ X4 + Y4

X1 ~~ Y1
X2 ~~ Y2
X3 ~~ Y3
X4 ~~ Y4
X5 ~~ Y5

FX ~~ FX
FY ~~ FY
FX ~~ FY

X1 ~~ X1
X2 ~~ X2
X3 ~~ X3
X4 ~~ X4
X5 ~~ X5
Y1 ~~ Y1
Y2 ~~ Y2
Y3 ~~ Y3
Y4 ~~ Y4
Y5 ~~ Y5

X1 ~ 1
X2 ~ 1
X3 ~ 1
X4 ~ 1
X5 ~ 1
Y1 ~ 1
Y2 ~ 1
Y3 ~ 1
Y4 ~ 1
Y5 ~ 1
"

model_dpm_free <- "
FX =~ NA*X2 + X3 + X4 + X5
FY =~ NA*Y2 + Y3 + Y4 + Y5
FX ~~ 1*FX
FY ~~ 1*FY
FX ~~ FY

FX ~~ X1 + Y1
FY ~~ X1 + Y1

X2 + Y2 ~ X1 + Y1
X3 + Y3 ~ X2 + Y2
X4 + Y4 ~ X3 + Y3
X5 + Y5 ~ X4 + Y4

X1 ~~ Y1
X2 ~~ Y2
X3 ~~ Y3
X4 ~~ Y4
X5 ~~ Y5

X1 ~~ X1
X2 ~~ X2
X3 ~~ X3
X4 ~~ X4
X5 ~~ X5
Y1 ~~ Y1
Y2 ~~ Y2
Y3 ~~ Y3
Y4 ~~ Y4
Y5 ~~ Y5

X1 ~ 1
X2 ~ 1
X3 ~ 1
X4 ~ 1
X5 ~ 1
Y1 ~ 1
Y2 ~ 1
Y3 ~ 1
Y4 ~ 1
Y5 ~ 1
"

model_riclpm <- "
RIX =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5
RIY =~ 1*Y1 + 1*Y2 + 1*Y3 + 1*Y4 + 1*Y5
RIX ~~ RIX
RIY ~~ RIY
RIX ~~ RIY

X1 ~~ 0*X1; X2 ~~ 0*X2; X3 ~~ 0*X3; X4 ~~ 0*X4; X5 ~~ 0*X5
Y1 ~~ 0*Y1; Y2 ~~ 0*Y2; Y3 ~~ 0*Y3; Y4 ~~ 0*Y4; Y5 ~~ 0*Y5

wX1 =~ 1*X1; wX2 =~ 1*X2; wX3 =~ 1*X3; wX4 =~ 1*X4; wX5 =~ 1*X5
wY1 =~ 1*Y1; wY2 =~ 1*Y2; wY3 =~ 1*Y3; wY4 =~ 1*Y4; wY5 =~ 1*Y5

RIX ~~ 0*wX1 + 0*wX2 + 0*wX3 + 0*wX4 + 0*wX5
RIX ~~ 0*wY1 + 0*wY2 + 0*wY3 + 0*wY4 + 0*wY5
RIY ~~ 0*wX1 + 0*wX2 + 0*wX3 + 0*wX4 + 0*wX5
RIY ~~ 0*wY1 + 0*wY2 + 0*wY3 + 0*wY4 + 0*wY5

wX1 ~~ wX1; wX2 ~~ wX2; wX3 ~~ wX3; wX4 ~~ wX4; wX5 ~~ wX5
wY1 ~~ wY1; wY2 ~~ wY2; wY3 ~~ wY3; wY4 ~~ wY4; wY5 ~~ wY5
wY1 ~~ wX1; wY2 ~~ wX2; wY3 ~~ wX3; wY4 ~~ wX4; wY5 ~~ wX5

wX2 ~ wX1 + wY1
wY2 ~ wX1 + wY1
wX3 ~ wX2 + wY2
wY3 ~ wX2 + wY2
wX4 ~ wX3 + wY3
wY4 ~ wX3 + wY3
wX5 ~ wX4 + wY4
wY5 ~ wX4 + wY4

X1 + X2 + X3 + X4 + X5 ~ mX*1
Y1 + Y2 + Y3 + Y4 + Y5 ~ mY*1
"

model_riclpm_free <- "
RIX =~ NA*X1 + X2 + X3 + X4 + X5
RIY =~ NA*Y1 + Y2 + Y3 + Y4 + Y5
RIX ~~ 1*RIX
RIY ~~ 1*RIY
RIX ~~ RIY

X1 ~~ 0*X1; X2 ~~ 0*X2; X3 ~~ 0*X3; X4 ~~ 0*X4; X5 ~~ 0*X5
Y1 ~~ 0*Y1; Y2 ~~ 0*Y2; Y3 ~~ 0*Y3; Y4 ~~ 0*Y4; Y5 ~~ 0*Y5

wX1 =~ 1*X1; wX2 =~ 1*X2; wX3 =~ 1*X3; wX4 =~ 1*X4; wX5 =~ 1*X5
wY1 =~ 1*Y1; wY2 =~ 1*Y2; wY3 =~ 1*Y3; wY4 =~ 1*Y4; wY5 =~ 1*Y5

RIX ~~ 0*wX1 + 0*wX2 + 0*wX3 + 0*wX4 + 0*wX5
RIX ~~ 0*wY1 + 0*wY2 + 0*wY3 + 0*wY4 + 0*wY5
RIY ~~ 0*wX1 + 0*wX2 + 0*wX3 + 0*wX4 + 0*wX5
RIY ~~ 0*wY1 + 0*wY2 + 0*wY3 + 0*wY4 + 0*wY5

wX1 ~~ wX1; wX2 ~~ wX2; wX3 ~~ wX3; wX4 ~~ wX4; wX5 ~~ wX5
wY1 ~~ wY1; wY2 ~~ wY2; wY3 ~~ wY3; wY4 ~~ wY4; wY5 ~~ wY5
wY1 ~~ wX1; wY2 ~~ wX2; wY3 ~~ wX3; wY4 ~~ wX4; wY5 ~~ wX5

wX2 ~ wX1 + wY1
wY2 ~ wX1 + wY1
wX3 ~ wX2 + wY2
wY3 ~ wX2 + wY2
wX4 ~ wX3 + wY3
wY4 ~ wX3 + wY3
wX5 ~ wX4 + wY4
wY5 ~ wX4 + wY4

X1 ~ 1; X2 ~ 1; X3 ~ 1; X4 ~ 1; X5 ~ 1
Y1 ~ 1; Y2 ~ 1; Y3 ~ 1; Y4 ~ 1; Y5 ~ 1
"

# =============================================================================
# Fitting helper
# =============================================================================
fit_model <- function(model, data = NULL, Sigma_full = NULL, Nobs = 500) {
  if (!is.null(data)) {
    lavaan::sem(model = model, data = data, estimator = "MLR")
  } else if (!is.null(Sigma_full)) {
    lavaan::sem(model = model, sample.cov = Sigma_full,
                sample.cov.rescale = FALSE, sample.nobs = Nobs)
  } else {
    stop("Either data or Sigma_full must be provided!")
  }
}

# =============================================================================
# Fit all models and return ONLY cross-lag effects; auto-add two covariate models
# =============================================================================
FitModels <- function(
    data = NULL,
    Sigma_full = NULL,
    Nobs = 500
){
  # --- 1) Detect covariates automatically from Sigma_full or data ------------
  all_names <- if (!is.null(Sigma_full)) {
    colnames(Sigma_full)
  } else if (!is.null(data)) {
    colnames(data)
  } else character(0)

  covariates <- all_names[grepl("C", all_names)]  # any name containing "C"

  # --- 2) Build the two covariate-augmented CLPMs (only if covariates exist) --
  if (length(covariates) > 0) {
    model_covariates <- paste0(
      model_clpm,
      "\nX1 + X2 + X3 + X4 + X5 ~ ",
      paste(covariates, collapse = " + "),
      "\nY1 + Y2 + Y3 + Y4 + Y5 ~ ",
      paste(covariates, collapse = " + ")
    )
    model_covariates_equal <- paste0(
      model_clpm,
      "\nX1 + X2 + X3 + X4 + X5 ~ ",
      paste(paste0("bX_", seq_along(covariates), "*", covariates), collapse = " + "),
      "\nY1 + Y2 + Y3 + Y4 + Y5 ~ ",
      paste(paste0("bY_", seq_along(covariates), "*", covariates), collapse = " + ")
    )
  } else {
    model_covariates <- model_clpm
    model_covariates_equal <- model_clpm
  }

  # --- 3) Fit all models ------------------------------------------------------
  fit_clpm              <- fit_model(model_clpm,              data, Sigma_full, Nobs)
  fit_clpm_lag2         <- fit_model(model_clpm_lag2,         data, Sigma_full, Nobs)
  fit_dpm               <- fit_model(model_dpm,               data, Sigma_full, Nobs)
  fit_dpm_free          <- fit_model(model_dpm_free,          data, Sigma_full, Nobs)
  fit_riclpm            <- fit_model(model_riclpm,            data, Sigma_full, Nobs)
  fit_riclpm_free       <- fit_model(model_riclpm_free,       data, Sigma_full, Nobs)
  fit_covariates        <- fit_model(model_covariates,        data, Sigma_full, Nobs)
  fit_covariates_equal  <- fit_model(model_covariates_equal,  data, Sigma_full, Nobs)

  # --- 4) Parameter tables ----------------------------------------------------
  PE_clpm        <- lavaan::parameterEstimates(fit_clpm)
  PE_clpm_lag2   <- lavaan::parameterEstimates(fit_clpm_lag2)
  PE_dpm         <- lavaan::parameterEstimates(fit_dpm)
  PE_dpm_free    <- lavaan::parameterEstimates(fit_dpm_free)
  PE_riclpm      <- lavaan::parameterEstimates(fit_riclpm)
  PE_riclpm_free <- lavaan::parameterEstimates(fit_riclpm_free)
  PE_cov         <- lavaan::parameterEstimates(fit_covariates)
  PE_cov_eq      <- lavaan::parameterEstimates(fit_covariates_equal)

  # --- 5) Selectors: keep ALL lagged (auto + cross); still lag-2 present for model_clpm_lag2 ----
  select_lagged_CLPM_like <- function(df) {
    df <- subset(df, op == "~" &
                   lhs %in% c(paste0("X", 2:5), paste0("Y", 2:5)) &
                   rhs %in% c(paste0("X", 1:4), paste0("Y", 1:4)))
    df  # keep BOTH auto (X~X, Y~Y) and cross (X~Y, Y~X)
  }
  select_lagged_RICLPM <- function(df) {
    tmp <- subset(df, op == "~" &
                    lhs %in% c(paste0("wX", 2:5), paste0("wY", 2:5)) &
                    rhs %in% c(paste0("wX", 1:4), paste0("wY", 1:4)))
    tmp$lhs <- sub("^wX","X", tmp$lhs); tmp$lhs <- sub("^wY","Y", tmp$lhs)
    tmp$rhs <- sub("^wX","X", tmp$rhs); tmp$rhs <- sub("^wY","Y", tmp$rhs)
    select_lagged_CLPM_like(tmp)
  }

  # --- 6) Extract and stack ---------------------------------------------------
  cl_tabs <- list(
    "model_clpm"               = select_lagged_CLPM_like(PE_clpm),
    "model_clpm_lag2"          = select_lagged_CLPM_like(PE_clpm_lag2),
    "model_dpm"                = select_lagged_CLPM_like(PE_dpm),
    "model_dpm_free"           = select_lagged_CLPM_like(PE_dpm_free),
    "model_riclpm"             = select_lagged_RICLPM(PE_riclpm),
    "model_riclpm_free"        = select_lagged_RICLPM(PE_riclpm_free),
    "model_covariates"         = select_lagged_CLPM_like(PE_cov),
    "model_covariates_equal"   = select_lagged_CLPM_like(PE_cov_eq)
  )

  keep_cols <- c("lhs","op","rhs","est","se","z")
  tag <- function(df, label) {
    out <- df[, intersect(colnames(df), keep_cols), drop = FALSE]
    for (nm in setdiff(keep_cols, names(out))) out[[nm]] <- NA_real_
    out$model <- label
    out
  }
  out <- do.call(rbind, Map(tag, cl_tabs, names(cl_tabs)))

  # --- 7) Annotate (now includes auto and cross) ------------------------------
  out$time <- as.numeric(sub("^[XY]([0-9]+)$", "\\1", out$rhs)) + 1L
  out$effect_type <- ifelse(substr(out$lhs,1,1) == substr(out$rhs,1,1),
                            "auto", "cross")
  out$effect <- ifelse(out$effect_type == "auto",
                       paste0(substr(out$rhs,1,1), "[t-1] -> ", substr(out$lhs,1,1), "[t]"),
                       paste0(substr(out$rhs,1,1), "[t-1] -> ", substr(out$lhs,1,1), "[t]"))

  out$model <- factor(out$model,
                      levels = c("model_covariates", "model_covariates_equal",
                                 "model_clpm","model_clpm_lag2",
                                 "model_dpm","model_dpm_free",
                                 "model_riclpm","model_riclpm_free"))
  out <- out[order(out$effect_type, out$time, out$model), ]

  # --- 8) Add Bias (relative to model_covariates; removes lag-2 automatically) ----
  temp <- out[out$model == "model_covariates", ]
  temp$estRef <- temp$est
  temp$est <- temp$se <- temp$z <- temp$model <- temp$time <- temp$effect <- temp$effect_type <- NULL
  out <- dplyr::left_join(out, temp)
  out <- out[!is.na(out$estRef), ]  # drops paths absent from reference (e.g., pure lag-2)
  rownames(out) <- NULL
  out$Bias <- out$est - out$estRef
  out$RelBiasPercent <- out$Bias / out$estRef * 100

  return(out)
}

get_coefs <- function(multi_list,
                      BetaBase){
  B_XY_List <- lapply(multi_list, FUN = function(M) M*BetaBase)
  return(B_XY_List)
}
