# This function calls the previous functions to run one replication of the simulation study.
#
# 1) It checks the inputs:
# - that B_scenarios is provided and contains all requested scenarios
# - that each scenario trajectory is a list of length T
# - that each B matrix is 2 x p, where p = k + number of interaction terms implied by k
# - that the A matrix is 2x2
# - that the Psi matrix is kxk (k = number of linear confounder parts)
#
# 2) It extracts the mean delta values (confounder effects) for X and Y at each time point for later output.
# 3) It simulates panel data using the provided parameters with simulate_panel_data_int().
# 4) It builds model strings for each model type.
# 5) It calls the appropriate function(s) to fit the requested models.
# 6) It calls the appropriate function(s) to extract gamma (cross-lagged), beta (autoregressive), and rho (residual correlation)
#    estimates, along with p-values and confidence intervals.
# 7) It stores the results in a list (one element per scenario) and returns a single bound data frame.
# ------------------------------------------------------------------------------------------------------------

if (!exists("n_interactions_from_k")) {
  n_interactions_from_k <- function(k) {
    if (is.na(k) || k < 2) return(0L)
    sum(vapply(2:k, function(m) choose(k, m), numeric(1)))
  }
}

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
  ci_level = 0.95,                                        # CI level for extracted parameters
  xgb_tuned = NULL,
  xgb_fit_profile = c("fast", "balanced", "thorough"),
  xgb_fit_overrides = NULL
){

  # match xgb fit profile argument
  xgb_fit_profile <- match.arg(xgb_fit_profile)

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

  # check A is 2 x 2
  if (!is.matrix(A) || any(dim(A) != c(2, 2)))
    stop("A must be a 2x2 matrix.")

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

    # extract mean deltas (confounder effects) for X and Y at each wave
    delta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    delta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    delta_vec   <- delta_X_vec

    # simulate panel data
    df <- tryCatch(

      # simulate panel data
      simulate_panel_data_int(
        N         = N,
        T         = T,
        A         = A,
        D_list    = B_list,
        Psi       = Psi,
        rho_extra = rho_extra
      ),
      error = function(e) NULL
    )

    # if simulation fails, return NA rows with the correct structure
    if (is.null(df)) {

      out_list[[j]] <- data.frame(

        # run info
        run      = rep(rep_id, T),
        occasion = 1:T,
        scenario = scen,

        # deltas = confounder effects
        delta    = delta_vec,
        delta_X  = delta_X_vec,
        delta_Y  = delta_Y_vec,

        # gammaXY = cross-lagged effect of X on Y
        estGammaXY_CLPM      = NA, pGammaXY_CLPM      = NA, ciL_GammaXY_CLPM      = NA, ciU_GammaXY_CLPM      = NA,
        estGammaXY_RI_CLPM   = NA, pGammaXY_RI_CLPM   = NA, ciL_GammaXY_RI_CLPM   = NA, ciU_GammaXY_RI_CLPM   = NA,
        estGammaXY_DPM       = NA, pGammaXY_DPM       = NA, ciL_GammaXY_DPM       = NA, ciU_GammaXY_DPM       = NA,
        estGammaXY_CLPM_Adj  = NA, pGammaXY_CLPM_Adj  = NA, ciL_GammaXY_CLPM_Adj  = NA, ciU_GammaXY_CLPM_Adj  = NA,
        estGammaXY_CLPM_LBCA = NA, pGammaXY_CLPM_LBCA = NA, ciL_GammaXY_CLPM_LBCA = NA, ciU_GammaXY_CLPM_LBCA = NA,
        estGammaXY_CLPM_XGB  = NA, pGammaXY_CLPM_XGB  = NA, ciL_GammaXY_CLPM_XGB  = NA, ciU_GammaXY_CLPM_XGB  = NA,

        # gammaYX = cross-lagged effect of Y on X
        estGammaYX_CLPM      = NA, pGammaYX_CLPM      = NA, ciL_GammaYX_CLPM      = NA, ciU_GammaYX_CLPM      = NA,
        estGammaYX_RI_CLPM   = NA, pGammaYX_RI_CLPM   = NA, ciL_GammaYX_RI_CLPM   = NA, ciU_GammaYX_RI_CLPM   = NA,
        estGammaYX_DPM       = NA, pGammaYX_DPM       = NA, ciL_GammaYX_DPM       = NA, ciU_GammaYX_DPM       = NA,
        estGammaYX_CLPM_Adj  = NA, pGammaYX_CLPM_Adj  = NA, ciL_GammaYX_CLPM_Adj  = NA, ciU_GammaYX_CLPM_Adj  = NA,
        estGammaYX_CLPM_LBCA = NA, pGammaYX_CLPM_LBCA = NA, ciL_GammaYX_CLPM_LBCA = NA, ciU_GammaYX_CLPM_LBCA = NA,
        estGammaYX_CLPM_XGB  = NA, pGammaYX_CLPM_XGB  = NA, ciL_GammaYX_CLPM_XGB  = NA, ciU_GammaYX_CLPM_XGB  = NA,

        # betaX = auto-regressive effect of X
        estBetaX_CLPM      = NA, pBetaX_CLPM      = NA, ciL_BetaX_CLPM      = NA, ciU_BetaX_CLPM      = NA,
        estBetaX_RI_CLPM   = NA, pBetaX_RI_CLPM   = NA, ciL_BetaX_RI_CLPM   = NA, ciU_BetaX_RI_CLPM   = NA,
        estBetaX_DPM       = NA, pBetaX_DPM       = NA, ciL_BetaX_DPM       = NA, ciU_BetaX_DPM       = NA,
        estBetaX_CLPM_Adj  = NA, pBetaX_CLPM_Adj  = NA, ciL_BetaX_CLPM_Adj  = NA, ciU_BetaX_CLPM_Adj  = NA,
        estBetaX_CLPM_LBCA = NA, pBetaX_CLPM_LBCA = NA, ciL_BetaX_CLPM_LBCA = NA, ciU_BetaX_CLPM_LBCA = NA,
        estBetaX_CLPM_XGB  = NA, pBetaX_CLPM_XGB  = NA, ciL_BetaX_CLPM_XGB  = NA, ciU_BetaX_CLPM_XGB  = NA,

        # betaY = auto-regressive effect of Y
        estBetaY_CLPM      = NA, pBetaY_CLPM      = NA, ciL_BetaY_CLPM      = NA, ciU_BetaY_CLPM      = NA,
        estBetaY_RI_CLPM   = NA, pBetaY_RI_CLPM   = NA, ciL_BetaY_RI_CLPM   = NA, ciU_BetaY_RI_CLPM   = NA,
        estBetaY_DPM       = NA, pBetaY_DPM       = NA, ciL_BetaY_DPM       = NA, ciU_BetaY_DPM       = NA,
        estBetaY_CLPM_Adj  = NA, pBetaY_CLPM_Adj  = NA, ciL_BetaY_CLPM_Adj  = NA, ciU_BetaY_CLPM_Adj  = NA,
        estBetaY_CLPM_LBCA = NA, pBetaY_CLPM_LBCA = NA, ciL_BetaY_CLPM_LBCA = NA, ciU_BetaY_CLPM_LBCA = NA,
        estBetaY_CLPM_XGB  = NA, pBetaY_CLPM_XGB  = NA, ciL_BetaY_CLPM_XGB  = NA, ciU_BetaY_CLPM_XGB  = NA,

        # rho = residual correlation
        estRho_CLPM      = NA, pRho_CLPM      = NA, ciL_Rho_CLPM      = NA, ciU_Rho_CLPM      = NA,
        estRho_RI_CLPM   = NA, pRho_RI_CLPM   = NA, ciL_Rho_RI_CLPM   = NA, ciU_Rho_RI_CLPM   = NA,
        estRho_DPM       = NA, pRho_DPM       = NA, ciL_Rho_DPM       = NA, ciU_Rho_DPM       = NA,
        estRho_CLPM_Adj  = NA, pRho_CLPM_Adj  = NA, ciL_Rho_CLPM_Adj  = NA, ciU_Rho_CLPM_Adj  = NA,
        estRho_CLPM_LBCA = NA, pRho_CLPM_LBCA = NA, ciL_Rho_CLPM_LBCA = NA, ciU_Rho_CLPM_LBCA = NA,
        estRho_CLPM_XGB  = NA, pRho_CLPM_XGB  = NA, ciL_Rho_CLPM_XGB  = NA, ciU_Rho_CLPM_XGB  = NA,

        # fail flags if model failed
        fail_CLPM      = TRUE,
        fail_RI_CLPM   = TRUE,
        fail_DPM       = TRUE,
        fail_CLPM_Adj  = TRUE,
        fail_CLPM_LBCA = TRUE,
        fail_CLPM_XGB  = TRUE,

        # error messages since sim failed
        err_CLPM      = "sim failed",
        err_RI_CLPM   = "sim failed",
        err_DPM       = "sim failed",
        err_CLPM_Adj  = "sim failed",
        err_CLPM_LBCA = "sim failed",
        err_CLPM_XGB  = "sim failed",
 
        # NA run marker
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
    res_xgb  <- if ("xgb"    %in% models_to_run) safe_fit_clpm_xgb(
      model_string = model_clpm,
      data = df,
      k = k,
      xgb_tuned = xgb_tuned,
      xgb_fit_profile = xgb_fit_profile,
      xgb_fit_overrides = xgb_fit_overrides
    ) else list(fit = NULL, err = NA_character_)

    # pull fitted objects
    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit
    fit_xgb      <- res_xgb$fit

    # extract lagged parameters
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",    ci_level = ci_level)
    lag_ric  <- extract_lagged_parameters(fit_ric,      T, "riclpm",  ci_level = ci_level)
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,     T, "dpm",     ci_level = ci_level)
    lag_adj  <- extract_lagged_parameters(fit_adj,      T, "clpm",    ci_level = ci_level)
    lag_lbca <- extract_lagged_parameters(fit_lbca,     T, "clpm",    ci_level = ci_level)
    lag_xgb  <- extract_lagged_parameters(fit_xgb,      T, "clpm",    ci_level = ci_level)

    # extract residual correlations
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm",   ci_level = ci_level)
    rho_ric  <- extract_rho_vec(fit_ric,      T, "riclpm", ci_level = ci_level)
    rho_dpm  <- extract_rho_vec(fit_dpm0,     T, "dpm",    ci_level = ci_level)
    rho_adj  <- extract_rho_vec(fit_adj,      T, "clpm",   ci_level = ci_level)
    rho_lbca <- extract_rho_vec(fit_lbca,     T, "clpm",   ci_level = ci_level)
    rho_xgb  <- extract_rho_vec(fit_xgb,      T, "clpm",   ci_level = ci_level)

    # assemble output rows
    out_list[[j]] <- data.frame(

      # run info
      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      # deltas = confounder effects
      delta    = delta_vec,
      delta_X  = delta_X_vec,
      delta_Y  = delta_Y_vec,

      # gammaXY = cross-lagged effect of X on Y
      estGammaXY_CLPM      = c(NA, lag_raw$xy$est),
      pGammaXY_CLPM        = c(NA, lag_raw$xy$p),
      ciL_GammaXY_CLPM     = c(NA, lag_raw$xy$ci.lower),
      ciU_GammaXY_CLPM     = c(NA, lag_raw$xy$ci.upper),

      estGammaXY_RI_CLPM   = c(NA, lag_ric$xy$est),
      pGammaXY_RI_CLPM     = c(NA, lag_ric$xy$p),
      ciL_GammaXY_RI_CLPM  = c(NA, lag_ric$xy$ci.lower),
      ciU_GammaXY_RI_CLPM  = c(NA, lag_ric$xy$ci.upper),

      estGammaXY_DPM       = c(NA, lag_dpm0$xy$est),
      pGammaXY_DPM         = c(NA, lag_dpm0$xy$p),
      ciL_GammaXY_DPM      = c(NA, lag_dpm0$xy$ci.lower),
      ciU_GammaXY_DPM      = c(NA, lag_dpm0$xy$ci.upper),

      estGammaXY_CLPM_Adj  = c(NA, lag_adj$xy$est),
      pGammaXY_CLPM_Adj    = c(NA, lag_adj$xy$p),
      ciL_GammaXY_CLPM_Adj = c(NA, lag_adj$xy$ci.lower),
      ciU_GammaXY_CLPM_Adj = c(NA, lag_adj$xy$ci.upper),

      estGammaXY_CLPM_LBCA  = c(NA, lag_lbca$xy$est),
      pGammaXY_CLPM_LBCA    = c(NA, lag_lbca$xy$p),
      ciL_GammaXY_CLPM_LBCA = c(NA, lag_lbca$xy$ci.lower),
      ciU_GammaXY_CLPM_LBCA = c(NA, lag_lbca$xy$ci.upper),

      estGammaXY_CLPM_XGB  = c(NA, lag_xgb$xy$est),
      pGammaXY_CLPM_XGB    = c(NA, lag_xgb$xy$p),
      ciL_GammaXY_CLPM_XGB = c(NA, lag_xgb$xy$ci.lower),
      ciU_GammaXY_CLPM_XGB = c(NA, lag_xgb$xy$ci.upper),

      # gamma YX = cross-lagged effect of Y on X
      estGammaYX_CLPM      = c(NA, lag_raw$yx$est),
      pGammaYX_CLPM        = c(NA, lag_raw$yx$p),
      ciL_GammaYX_CLPM     = c(NA, lag_raw$yx$ci.lower),
      ciU_GammaYX_CLPM     = c(NA, lag_raw$yx$ci.upper),

      estGammaYX_RI_CLPM   = c(NA, lag_ric$yx$est),
      pGammaYX_RI_CLPM     = c(NA, lag_ric$yx$p),
      ciL_GammaYX_RI_CLPM  = c(NA, lag_ric$yx$ci.lower),
      ciU_GammaYX_RI_CLPM  = c(NA, lag_ric$yx$ci.upper),

      estGammaYX_DPM       = c(NA, lag_dpm0$yx$est),
      pGammaYX_DPM         = c(NA, lag_dpm0$yx$p),
      ciL_GammaYX_DPM      = c(NA, lag_dpm0$yx$ci.lower),
      ciU_GammaYX_DPM      = c(NA, lag_dpm0$yx$ci.upper),

      estGammaYX_CLPM_Adj  = c(NA, lag_adj$yx$est),
      pGammaYX_CLPM_Adj    = c(NA, lag_adj$yx$p),
      ciL_GammaYX_CLPM_Adj = c(NA, lag_adj$yx$ci.lower),
      ciU_GammaYX_CLPM_Adj = c(NA, lag_adj$yx$ci.upper),

      estGammaYX_CLPM_LBCA  = c(NA, lag_lbca$yx$est),
      pGammaYX_CLPM_LBCA    = c(NA, lag_lbca$yx$p),
      ciL_GammaYX_CLPM_LBCA = c(NA, lag_lbca$yx$ci.lower),
      ciU_GammaYX_CLPM_LBCA = c(NA, lag_lbca$yx$ci.upper),

      estGammaYX_CLPM_XGB  = c(NA, lag_xgb$yx$est),
      pGammaYX_CLPM_XGB    = c(NA, lag_xgb$yx$p),
      ciL_GammaYX_CLPM_XGB = c(NA, lag_xgb$yx$ci.lower),
      ciU_GammaYX_CLPM_XGB = c(NA, lag_xgb$yx$ci.upper),

      # betaX = auto-regressive effect of X
      estBetaX_CLPM      = c(NA, lag_raw$ar_x$est),
      pBetaX_CLPM        = c(NA, lag_raw$ar_x$p),
      ciL_BetaX_CLPM     = c(NA, lag_raw$ar_x$ci.lower),
      ciU_BetaX_CLPM     = c(NA, lag_raw$ar_x$ci.upper),

      estBetaX_RI_CLPM   = c(NA, lag_ric$ar_x$est),
      pBetaX_RI_CLPM     = c(NA, lag_ric$ar_x$p),
      ciL_BetaX_RI_CLPM  = c(NA, lag_ric$ar_x$ci.lower),
      ciU_BetaX_RI_CLPM  = c(NA, lag_ric$ar_x$ci.upper),

      estBetaX_DPM       = c(NA, lag_dpm0$ar_x$est),
      pBetaX_DPM         = c(NA, lag_dpm0$ar_x$p),
      ciL_BetaX_DPM      = c(NA, lag_dpm0$ar_x$ci.lower),
      ciU_BetaX_DPM      = c(NA, lag_dpm0$ar_x$ci.upper),

      estBetaX_CLPM_Adj  = c(NA, lag_adj$ar_x$est),
      pBetaX_CLPM_Adj    = c(NA, lag_adj$ar_x$p),
      ciL_BetaX_CLPM_Adj = c(NA, lag_adj$ar_x$ci.lower),
      ciU_BetaX_CLPM_Adj = c(NA, lag_adj$ar_x$ci.upper),

      estBetaX_CLPM_LBCA  = c(NA, lag_lbca$ar_x$est),
      pBetaX_CLPM_LBCA    = c(NA, lag_lbca$ar_x$p),
      ciL_BetaX_CLPM_LBCA = c(NA, lag_lbca$ar_x$ci.lower),
      ciU_BetaX_CLPM_LBCA = c(NA, lag_lbca$ar_x$ci.upper),

      estBetaX_CLPM_XGB  = c(NA, lag_xgb$ar_x$est),
      pBetaX_CLPM_XGB    = c(NA, lag_xgb$ar_x$p),
      ciL_BetaX_CLPM_XGB = c(NA, lag_xgb$ar_x$ci.lower),
      ciU_BetaX_CLPM_XGB = c(NA, lag_xgb$ar_x$ci.upper),

      # betaY = auto-regressive effect of Y
      estBetaY_CLPM      = c(NA, lag_raw$ar_y$est),
      pBetaY_CLPM        = c(NA, lag_raw$ar_y$p),
      ciL_BetaY_CLPM     = c(NA, lag_raw$ar_y$ci.lower),
      ciU_BetaY_CLPM     = c(NA, lag_raw$ar_y$ci.upper),

      estBetaY_RI_CLPM   = c(NA, lag_ric$ar_y$est),
      pBetaY_RI_CLPM     = c(NA, lag_ric$ar_y$p),
      ciL_BetaY_RI_CLPM  = c(NA, lag_ric$ar_y$ci.lower),
      ciU_BetaY_RI_CLPM  = c(NA, lag_ric$ar_y$ci.upper),

      estBetaY_DPM       = c(NA, lag_dpm0$ar_y$est),
      pBetaY_DPM         = c(NA, lag_dpm0$ar_y$p),
      ciL_BetaY_DPM      = c(NA, lag_dpm0$ar_y$ci.lower),
      ciU_BetaY_DPM      = c(NA, lag_dpm0$ar_y$ci.upper),

      estBetaY_CLPM_Adj  = c(NA, lag_adj$ar_y$est),
      pBetaY_CLPM_Adj    = c(NA, lag_adj$ar_y$p),
      ciL_BetaY_CLPM_Adj = c(NA, lag_adj$ar_y$ci.lower),
      ciU_BetaY_CLPM_Adj = c(NA, lag_adj$ar_y$ci.upper),

      estBetaY_CLPM_LBCA  = c(NA, lag_lbca$ar_y$est),
      pBetaY_CLPM_LBCA    = c(NA, lag_lbca$ar_y$p),
      ciL_BetaY_CLPM_LBCA = c(NA, lag_lbca$ar_y$ci.lower),
      ciU_BetaY_CLPM_LBCA = c(NA, lag_lbca$ar_y$ci.upper),

      estBetaY_CLPM_XGB  = c(NA, lag_xgb$ar_y$est),
      pBetaY_CLPM_XGB    = c(NA, lag_xgb$ar_y$p),
      ciL_BetaY_CLPM_XGB = c(NA, lag_xgb$ar_y$ci.lower),
      ciU_BetaY_CLPM_XGB = c(NA, lag_xgb$ar_y$ci.upper),

      # rho = residual correlation between X and Y
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

      # fit failures
      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),
      fail_CLPM_XGB  = is.null(fit_xgb),

      # error messages
      err_CLPM      = rep(res_clpm$err,  T),
      err_RI_CLPM   = rep(res_ric$err,   T),
      err_DPM       = rep(res_dpm0$err,  T),
      err_CLPM_Adj  = rep(res_adj$err,   T),
      err_CLPM_LBCA = rep(res_lbca$err,  T),
      err_CLPM_XGB  = rep(res_xgb$err,   T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy$est, lag_ric$xy$est, lag_dpm0$xy$est,
        lag_adj$xy$est, lag_lbca$xy$est, lag_xgb$xy$est
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}
