# This script contains a function that runs all other functions for one replication of the study
# this function should contain all the arguments that are defined in the previous sections
# ------------------------------------------------------------------------------------------------

run_one_rep_study <- function(
  rep_id,                                                 # replication index (set by outer loop)
  N,                                                      # sample size
  T,                                                      # number of waves
  k,                                                      # number of linear confounder parts
  scenarios,                                              # scenario names
  B_scenarios,                                            # named list: B_scenarios[[scen]] is a list of length T of 2 x p matrices
  A,                                                      # 2×2 autoregressive + cross-lag matrix
  Psi,                                                    # k×k confounder covariance (linear parts only)
  rho_extra,                                              # extra covariance added to X,Y each wave
  models_to_run,                                          # e.g. c("clpm","riclpm","dpm","lbca","adj","xgb")
  base_seed = 1234,                                       # base seed
  ci_level = 0.95                                         # CI level for extracted parameters
){

  # set seed for this replication
  set.seed(base_seed + rep_id)

  # expected number of interaction terms from k linear confounders
  k_non <- n_interactions_from_k(k)

  # expected total number of confounder columns (linear + interactions)
  p_exp <- k + k_non

  # check beta scenarios were provided
  if (is.null(B_scenarios))
    stop("B_scenarios is NULL. Please provide pre-defined beta trajectories via B_scenarios.")

  # check that scenarios exist in B_scenarios
  missing_scens <- setdiff(scenarios, names(B_scenarios))
  if (length(missing_scens) > 0)
    stop("Scenario(s) not found in B_scenarios: ", paste(missing_scens, collapse = ", "))

  # check Psi is k x k
  if (!is.matrix(Psi) || any(dim(Psi) != c(k, k)))
    stop("Mismatch: Psi must be k x k with k = ", k, ".")

  # prepare output list for each scenario
  out_list <- vector("list", length(scenarios))

  # loop over the scenarios
  for (j in seq_along(scenarios)) {

    # scenario name
    scen <- scenarios[j]

    # take the beta trajectory for this scenario
    B_list <- B_scenarios[[scen]]

    # check B_list is a list
    if (!is.list(B_list))
      stop("B_scenarios[['", scen, "']] must be a list of length T.")

    # check B_list has length T
    if (length(B_list) != T)
      stop("B_scenarios[['", scen, "']] has length ", length(B_list), " but T = ", T, ".")

    # check each matrix is 2 x p_exp
    for (t in 1:T) {

      # check element is a matrix
      if (!is.matrix(B_list[[t]]))
        stop("B_scenarios[['", scen, "']][[", t, "]] is not a matrix.")

      # check rows are X and Y
      if (nrow(B_list[[t]]) != 2L)
        stop("B_scenarios[['", scen, "']][[", t, "]] must have 2 rows (X,Y).")

      # check total number of confounder columns
      if (ncol(B_list[[t]]) != p_exp)
        stop("B_scenarios[['", scen, "']][[", t, "]] has p = ", ncol(B_list[[t]]),
             " columns, but expected p = k + #interactions = ", p_exp, ".")
    }

    # extract mean betas
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec

    # simulate panel data
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

    # if simulation fails, return NA rows with the correct structure
    if (is.null(df)) {

      out_list[[j]] <- data.frame(
        run      = rep(rep_id, T),
        occasion = 1:T,
        scenario = scen,
        beta     = beta_vec,
        beta_X   = beta_X_vec,
        beta_Y   = beta_Y_vec,

        # XY
        estXY_CLPM      = NA, pXY_CLPM      = NA, ciL_XY_CLPM      = NA, ciU_XY_CLPM      = NA,
        estXY_RI_CLPM   = NA, pXY_RI_CLPM   = NA, ciL_XY_RI_CLPM   = NA, ciU_XY_RI_CLPM   = NA,
        estXY_DPM       = NA, pXY_DPM       = NA, ciL_XY_DPM       = NA, ciU_XY_DPM       = NA,
        estXY_CLPM_Adj  = NA, pXY_CLPM_Adj  = NA, ciL_XY_CLPM_Adj  = NA, ciU_XY_CLPM_Adj  = NA,
        estXY_CLPM_LBCA = NA, pXY_CLPM_LBCA = NA, ciL_XY_CLPM_LBCA = NA, ciU_XY_CLPM_LBCA = NA,
        estXY_CLPM_XGB  = NA, pXY_CLPM_XGB  = NA, ciL_XY_CLPM_XGB  = NA, ciU_XY_CLPM_XGB  = NA,

        # YX
        estYX_CLPM      = NA, pYX_CLPM      = NA, ciL_YX_CLPM      = NA, ciU_YX_CLPM      = NA,
        estYX_RI_CLPM   = NA, pYX_RI_CLPM   = NA, ciL_YX_RI_CLPM   = NA, ciU_YX_RI_CLPM   = NA,
        estYX_DPM       = NA, pYX_DPM       = NA, ciL_YX_DPM       = NA, ciU_YX_DPM       = NA,
        estYX_CLPM_Adj  = NA, pYX_CLPM_Adj  = NA, ciL_YX_CLPM_Adj  = NA, ciU_YX_CLPM_Adj  = NA,
        estYX_CLPM_LBCA = NA, pYX_CLPM_LBCA = NA, ciL_YX_CLPM_LBCA = NA, ciU_YX_CLPM_LBCA = NA,
        estYX_CLPM_XGB  = NA, pYX_CLPM_XGB  = NA, ciL_YX_CLPM_XGB  = NA, ciU_YX_CLPM_XGB  = NA,

        # AX
        estA_CLPM      = NA, pA_CLPM      = NA, ciL_A_CLPM      = NA, ciU_A_CLPM      = NA,
        estA_RI_CLPM   = NA, pA_RI_CLPM   = NA, ciL_A_RI_CLPM   = NA, ciU_A_RI_CLPM   = NA,
        estA_DPM       = NA, pA_DPM       = NA, ciL_A_DPM       = NA, ciU_A_DPM       = NA,
        estA_CLPM_Adj  = NA, pA_CLPM_Adj  = NA, ciL_A_CLPM_Adj  = NA, ciU_A_CLPM_Adj  = NA,
        estA_CLPM_LBCA = NA, pA_CLPM_LBCA = NA, ciL_A_CLPM_LBCA = NA, ciU_A_CLPM_LBCA = NA,
        estA_CLPM_XGB  = NA, pA_CLPM_XGB  = NA, ciL_A_CLPM_XGB  = NA, ciU_A_CLPM_XGB  = NA,

        # AY
        estAY_CLPM      = NA, pAY_CLPM      = NA, ciL_AY_CLPM      = NA, ciU_AY_CLPM      = NA,
        estAY_RI_CLPM   = NA, pAY_RI_CLPM   = NA, ciL_AY_RI_CLPM   = NA, ciU_AY_RI_CLPM   = NA,
        estAY_DPM       = NA, pAY_DPM       = NA, ciL_AY_DPM       = NA, ciU_AY_DPM       = NA,
        estAY_CLPM_Adj  = NA, pAY_CLPM_Adj  = NA, ciL_AY_CLPM_Adj  = NA, ciU_AY_CLPM_Adj  = NA,
        estAY_CLPM_LBCA = NA, pAY_CLPM_LBCA = NA, ciL_AY_CLPM_LBCA = NA, ciU_AY_CLPM_LBCA = NA,
        estAY_CLPM_XGB  = NA, pAY_CLPM_XGB  = NA, ciL_AY_CLPM_XGB  = NA, ciU_AY_CLPM_XGB  = NA,

        # Rho
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
    res_clpm <- if ("clpm"   %in% models_to_run) safe_fit_clpm(model_clpm, df)           else list(fit = NULL, err = NA_character_)
    res_ric  <- if ("riclpm" %in% models_to_run) safe_fit_riclpm(model_riclpm, df)       else list(fit = NULL, err = NA_character_)
    res_dpm0 <- if ("dpm"    %in% models_to_run) safe_fit_dpm(model_dpm, df)             else list(fit = NULL, err = NA_character_)
    res_adj  <- if ("adj"    %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df) else list(fit = NULL, err = NA_character_)

    # fit LBCA and XGB using only the linear confounders (c1..ck)
    res_lbca <- if ("lbca"   %in% models_to_run) safe_fit_clpm_resid(model_clpm, df, k)  else list(fit = NULL, err = NA_character_)
    res_xgb  <- if ("xgb"    %in% models_to_run) safe_fit_clpm_xgb(model_clpm, df, k)    else list(fit = NULL, err = NA_character_)

    # pull fitted objects
    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit
    fit_xgb      <- res_xgb$fit

    # extract lagged parameters
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",    ci_level = ci_level)
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm",  ci_level = ci_level)
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm",     ci_level = ci_level)
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm",    ci_level = ci_level)
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm",    ci_level = ci_level)
    lag_xgb  <- extract_lagged_parameters(fit_xgb,       T, "clpm",    ci_level = ci_level)

    # extract residual correlations
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm",   ci_level = ci_level)
    rho_ric  <- extract_rho_vec(fit_ric,      T, "riclpm", ci_level = ci_level)
    rho_dpm  <- extract_rho_vec(fit_dpm0,     T, "dpm",    ci_level = ci_level)
    rho_adj  <- extract_rho_vec(fit_adj,      T, "clpm",   ci_level = ci_level)
    rho_lbca <- extract_rho_vec(fit_lbca,     T, "clpm",   ci_level = ci_level)
    rho_xgb  <- extract_rho_vec(fit_xgb,      T, "clpm",   ci_level = ci_level)

    # assemble output rows
    out_list[[j]] <- data.frame(
      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      beta     = beta_vec,
      beta_X   = beta_X_vec,
      beta_Y   = beta_Y_vec,

      # XY
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

      # YX
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

      # AX
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

      # AY
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

      # Rho
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

      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),
      fail_CLPM_XGB  = is.null(fit_xgb),

      err_CLPM      = rep(res_clpm$err,  T),
      err_RI_CLPM   = rep(res_ric$err,   T),
      err_DPM       = rep(res_dpm0$err,  T),
      err_CLPM_Adj  = rep(res_adj$err,   T),
      err_CLPM_LBCA = rep(res_lbca$err,  T),
      err_CLPM_XGB  = rep(res_xgb$err,   T),

      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy$est, lag_ric$xy$est, lag_dpm0$xy$est,
        lag_adj$xy$est, lag_lbca$xy$est, lag_xgb$xy$est
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}
