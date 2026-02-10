# This function calls the previous functions to run one replication of the simulation study. 
#
# 1) It checks the inputs:
# - that each delta matrix is 2xk
# - that each scenario exists in D_scenarios
# - that each delta trajectory has length T
# - that the A matrix is 2x2
# - that the Psi matrix is kxk
#
# 2) It extracts the mean delta values for X and Y at each time point for later output.
# 3) It simulates panel data using the provided parameters with simulate_panel_data().
# 4) It builds model strings for each model type.
# 5) It calls the appropriate function to fit the models.
# 6) It checks convergence and properness (improper solutions are flagged separately from non-convergence).
# 7) It calls the appropriate function to extract parameter estimates and confidence intervals.
# 8) It stores the results in a list.
# ------------------------------------------------------------------------------------------------------------

run_one_rep_study <- function(
  rep_id,                                                                               # replication index (set by outer loop)
  N,                                                                                    # sample size
  T,                                                                                    # number of waves
  k,                                                                                    # number of confounders
  scenarios,                                                                            # scenario names
  D_scenarios,                                                                          # named list of delta trajectories
  A,                                                                                    # 2×2 autoregressive (beta) + cross-lagged (gamma) matrix
  Psi,                                                                                  # k×k confounder covariance
  rho_extra = 0.15,                                                                     # extra covariance added to X,Y at time t
  models_to_run = c("clpm","adj","riclpm","dpm","riclpm_free","dpm_free","lbca","bca_riclpm","bca_dpm"),
  base_seed = 1234,                                                                     # base seed
  ci_level = 0.95                                                                       # CI level for extracted parameters
){

  # set seed for this replication
  set.seed(base_seed + rep_id)

  # check delta scenarios were provided
  if (is.null(D_scenarios))
    stop("D_scenarios is NULL. Please provide pre-defined delta trajectories via D_scenarios.")

  # check that scenarios exist in D_scenarios
  missing_scens <- setdiff(scenarios, names(D_scenarios))
  if (length(missing_scens) > 0)
    stop("Scenario(s) not found in D_scenarios: ", paste(missing_scens, collapse = ", "))

  # check A is 2 x 2
  if (!is.matrix(A) || any(dim(A) != c(2, 2)))
    stop("A must be a 2x2 matrix.")

  # check Psi is k x k
  if (!is.matrix(Psi) || any(dim(Psi) != c(k, k)))
    stop("Psi must be a k x k matrix with k = ", k, ".")

  # prepare output list for each scenario
  out_list <- vector("list", length(scenarios))

  # loop over the scenarios
  for (j in seq_along(scenarios)) {

    # scenario name
    scen <- scenarios[j]

    # take the delta trajectory for this scenario
    D_list <- D_scenarios[[scen]]

    # check the trajectory has length T
    if (!is.list(D_list) || length(D_list) != T)
      stop("D_scenarios[['", scen, "']] must be a list of length T.")

    # check each matrix is 2 x k
    bad_dims <- vapply(D_list, function(Dt) {
      !is.matrix(Dt) || any(dim(Dt) != c(2, k))
    }, logical(1))

    # stop if any bad dimensions
    if (any(bad_dims))
      stop("One or more D matrices in scenario '", scen, "' are not 2 x k with k = ", k, ".")

    # extract mean deltas
    delta_X_vec <- sapply(D_list, function(Dt) mean(Dt[1, ]))
    delta_Y_vec <- sapply(D_list, function(Dt) mean(Dt[2, ]))
    delta_vec   <- delta_X_vec

    # simulate panel data
    df <- tryCatch(
      simulate_panel_data(
        N         = N,
        T         = T,
        A         = A,
        D_list    = D_list,
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

        # true parameters (raw)
        betaX_true  = rep(A[1,1], T),
        betaY_true  = rep(A[2,2], T),
        gammaXY_true = rep(A[2,1], T),
        gammaYX_true = rep(A[1,2], T),

        # delta (confounder effects) 
        delta    = delta_vec,
        delta_X  = delta_X_vec,
        delta_Y  = delta_Y_vec,

        # cross-lagged effects = gamma
        estGammaXY_CLPM            = NA,
        estGammaXY_RI_CLPM         = NA,
        estGammaXY_DPM             = NA,
        estGammaXY_RI_CLPM_Free    = NA,
        estGammaXY_DPM_Free        = NA,
        estGammaXY_CLPM_Adj        = NA,
        estGammaXY_CLPM_LBCA       = NA,
        estGammaXY_RI_CLPM_BCA     = NA,
        estGammaXY_DPM_BCA         = NA,

        estGammaYX_CLPM            = NA,
        estGammaYX_RI_CLPM         = NA,
        estGammaYX_DPM             = NA,
        estGammaYX_RI_CLPM_Free    = NA,
        estGammaYX_DPM_Free        = NA,
        estGammaYX_CLPM_Adj        = NA,
        estGammaYX_CLPM_LBCA       = NA,
        estGammaYX_RI_CLPM_BCA     = NA,
        estGammaYX_DPM_BCA         = NA,

        # autoregressive effects = beta
        estBetaX_CLPM            = NA,
        estBetaX_RI_CLPM         = NA,
        estBetaX_DPM             = NA,
        estBetaX_RI_CLPM_Free    = NA,
        estBetaX_DPM_Free        = NA,
        estBetaX_CLPM_Adj        = NA,
        estBetaX_CLPM_LBCA       = NA,
        estBetaX_RI_CLPM_BCA     = NA,
        estBetaX_DPM_BCA         = NA,

        estBetaY_CLPM            = NA,
        estBetaY_RI_CLPM         = NA,
        estBetaY_DPM             = NA,
        estBetaY_RI_CLPM_Free    = NA,
        estBetaY_DPM_Free        = NA,
        estBetaY_CLPM_Adj        = NA,
        estBetaY_CLPM_LBCA       = NA,
        estBetaY_RI_CLPM_BCA     = NA,
        estBetaY_DPM_BCA         = NA,

        # residual correlations from CLPM = rho
        estRho_CLPM = NA,

        # model fit failures (non-convergence or fit object missing)
        fail_CLPM           = TRUE,
        fail_RI_CLPM        = TRUE,
        fail_DPM            = TRUE,
        fail_RI_CLPM_Free   = TRUE,
        fail_DPM_Free       = TRUE,
        fail_CLPM_Adj       = TRUE,
        fail_CLPM_LBCA      = TRUE,
        fail_RI_CLPM_BCA    = TRUE,
        fail_DPM_BCA        = TRUE,

        # improper flags (only relevant if converged)
        improper_CLPM           = NA,
        improper_RI_CLPM        = NA,
        improper_DPM            = NA,
        improper_RI_CLPM_Free   = NA,
        improper_DPM_Free       = NA,
        improper_CLPM_Adj       = NA,
        improper_CLPM_LBCA      = NA,
        improper_RI_CLPM_BCA    = NA,
        improper_DPM_BCA        = NA,

        # improper reasons (collapsed to a single string)
        improper_reason_CLPM           = NA,
        improper_reason_RI_CLPM        = NA,
        improper_reason_DPM            = NA,
        improper_reason_RI_CLPM_Free   = NA,
        improper_reason_DPM_Free       = NA,
        improper_reason_CLPM_Adj       = NA,
        improper_reason_CLPM_LBCA      = NA,
        improper_reason_RI_CLPM_BCA    = NA,
        improper_reason_DPM_BCA        = NA,

        # error messages if sim failed
        err_CLPM           = "sim failed",
        err_RI_CLPM        = "sim failed",
        err_DPM            = "sim failed",
        err_RI_CLPM_Free   = "sim failed",
        err_DPM_Free       = "sim failed",
        err_CLPM_Adj       = "sim failed",
        err_CLPM_LBCA      = "sim failed",
        err_RI_CLPM_BCA    = "sim failed",
        err_DPM_BCA        = "sim failed",

        # NA run marker
        is_na_run = 1L
      )

      next
    }

    # build model strings
    model_clpm             <- build_clpm(T)
    model_riclpm           <- build_riclpm(T)
    model_dpm              <- build_dpm(T)
    model_riclpm_free      <- build_riclpm_free_ri_loadings(T)
    model_dpm_free         <- build_dpm_free_loadings(T)
    model_clpm_with_Cs     <- build_clpm_with_Cs(T, k)

    # fit the models safely
    res_clpm      <- if ("clpm"        %in% models_to_run) safe_fit_clpm(model_clpm, df)                 else list(fit = NULL, err = NA_character_)
    res_adj       <- if ("adj"         %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df)       else list(fit = NULL, err = NA_character_)
    res_ric       <- if ("riclpm"      %in% models_to_run) safe_fit_riclpm(model_riclpm, df)             else list(fit = NULL, err = NA_character_)
    res_dpm0      <- if ("dpm"         %in% models_to_run) safe_fit_dpm(model_dpm, df)                   else list(fit = NULL, err = NA_character_)
    res_ric_free  <- if ("riclpm_free" %in% models_to_run) safe_fit_riclpm_free_loadings(model_riclpm_free, df) else list(fit = NULL, err = NA_character_)
    res_dpm_free  <- if ("dpm_free"    %in% models_to_run) safe_fit_dpm_free_loadings(model_dpm_free, df)       else list(fit = NULL, err = NA_character_)
    res_lbca      <- if ("lbca"        %in% models_to_run) safe_fit_clpm_resid(model_clpm, df)           else list(fit = NULL, err = NA_character_)
    res_bca_ric   <- if ("bca_riclpm"  %in% models_to_run) safe_fit_riclpm_resid(model_riclpm, df)       else list(fit = NULL, err = NA_character_)
    res_bca_dpm   <- if ("bca_dpm"     %in% models_to_run) safe_fit_dpm_resid(model_dpm, df)             else list(fit = NULL, err = NA_character_)

    # pull fitted objects
    fit_clpm_raw  <- res_clpm$fit
    fit_adj       <- res_adj$fit
    fit_ric       <- res_ric$fit
    fit_dpm0      <- res_dpm0$fit
    fit_ric_free  <- res_ric_free$fit
    fit_dpm_free  <- res_dpm_free$fit
    fit_lbca      <- res_lbca$fit
    fit_bca_ric   <- res_bca_ric$fit
    fit_bca_dpm   <- res_bca_dpm$fit

    # check convergence and properness
    chk_clpm     <- check_convergence_and_properness(fit_clpm_raw)
    chk_adj      <- check_convergence_and_properness(fit_adj)
    chk_ric      <- check_convergence_and_properness(fit_ric)
    chk_dpm0     <- check_convergence_and_properness(fit_dpm0)
    chk_ric_free <- check_convergence_and_properness(fit_ric_free)
    chk_dpm_free <- check_convergence_and_properness(fit_dpm_free)
    chk_lbca     <- check_convergence_and_properness(fit_lbca)
    chk_bca_ric  <- check_convergence_and_properness(fit_bca_ric)
    chk_bca_dpm  <- check_convergence_and_properness(fit_bca_dpm)

    # decide which fits are usable for extraction (must converge and be proper)
    fit_clpm_use     <- if (!is.null(fit_clpm_raw) && isTRUE(chk_clpm$converged)     && isTRUE(chk_clpm$proper))     fit_clpm_raw  else NULL
    fit_adj_use      <- if (!is.null(fit_adj)      && isTRUE(chk_adj$converged)      && isTRUE(chk_adj$proper))      fit_adj       else NULL
    fit_ric_use      <- if (!is.null(fit_ric)      && isTRUE(chk_ric$converged)      && isTRUE(chk_ric$proper))      fit_ric       else NULL
    fit_dpm0_use     <- if (!is.null(fit_dpm0)     && isTRUE(chk_dpm0$converged)     && isTRUE(chk_dpm0$proper))     fit_dpm0      else NULL
    fit_ric_free_use <- if (!is.null(fit_ric_free) && isTRUE(chk_ric_free$converged) && isTRUE(chk_ric_free$proper)) fit_ric_free  else NULL
    fit_dpm_free_use <- if (!is.null(fit_dpm_free) && isTRUE(chk_dpm_free$converged) && isTRUE(chk_dpm_free$proper)) fit_dpm_free  else NULL
    fit_lbca_use     <- if (!is.null(fit_lbca)     && isTRUE(chk_lbca$converged)     && isTRUE(chk_lbca$proper))     fit_lbca      else NULL
    fit_bca_ric_use  <- if (!is.null(fit_bca_ric)  && isTRUE(chk_bca_ric$converged)  && isTRUE(chk_bca_ric$proper))  fit_bca_ric   else NULL
    fit_bca_dpm_use  <- if (!is.null(fit_bca_dpm)  && isTRUE(chk_bca_dpm$converged)  && isTRUE(chk_bca_dpm$proper))  fit_bca_dpm   else NULL

    # extract lagged parameters
    lag_raw      <- extract_lagged_parameters(fit_clpm_use,     T, "clpm",   ci_level = ci_level)
    lag_adj      <- extract_lagged_parameters(fit_adj_use,      T, "clpm",   ci_level = ci_level)
    lag_ric      <- extract_lagged_parameters(fit_ric_use,      T, "riclpm", ci_level = ci_level)
    lag_dpm0     <- extract_lagged_parameters(fit_dpm0_use,     T, "dpm",    ci_level = ci_level)
    lag_ric_free <- extract_lagged_parameters(fit_ric_free_use, T, "riclpm", ci_level = ci_level)
    lag_dpm_free <- extract_lagged_parameters(fit_dpm_free_use, T, "dpm",    ci_level = ci_level)
    lag_lbca     <- extract_lagged_parameters(fit_lbca_use,     T, "clpm",   ci_level = ci_level)
    lag_bca_ric  <- extract_lagged_parameters(fit_bca_ric_use,  T, "riclpm", ci_level = ci_level)
    lag_bca_dpm  <- extract_lagged_parameters(fit_bca_dpm_use,  T, "dpm",    ci_level = ci_level)

    # extract residual correlations from CLPM (only if proper)
    rho_clpm <- extract_rho_vec(fit_clpm_use, T, "clpm")

    # helper to collapse reasons into one string (or NA)
    collapse_reasons <- function(chk) {
      if (is.null(chk$reasons) || length(chk$reasons) == 0) return(NA_character_)
      if (length(chk$reasons) == 1 && is.na(chk$reasons)) return(NA_character_)
      paste(chk$reasons, collapse = " | ")
    }

    # assemble output rows
    out_list[[j]] <- data.frame(

      # run info
      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      # true parameters (raw)
      betaX_true   = rep(A[1,1], T),
      betaY_true   = rep(A[2,2], T),
      gammaXY_true = rep(A[2,1], T),
      gammaYX_true = rep(A[1,2], T),

      # delta (confounder effects) 
      delta    = delta_vec,
      delta_X  = delta_X_vec,
      delta_Y  = delta_Y_vec,

      # cross-lagged effects = gamma
      estGammaXY_CLPM            = c(NA, lag_raw$xy[, "est"]),
      estGammaXY_RI_CLPM         = c(NA, lag_ric$xy[, "est"]),
      estGammaXY_DPM             = c(NA, lag_dpm0$xy[, "est"]),
      estGammaXY_RI_CLPM_Free    = c(NA, lag_ric_free$xy[, "est"]),
      estGammaXY_DPM_Free        = c(NA, lag_dpm_free$xy[, "est"]),
      estGammaXY_CLPM_Adj        = c(NA, lag_adj$xy[, "est"]),
      estGammaXY_CLPM_LBCA       = c(NA, lag_lbca$xy[, "est"]),
      estGammaXY_RI_CLPM_BCA     = c(NA, lag_bca_ric$xy[, "est"]),
      estGammaXY_DPM_BCA         = c(NA, lag_bca_dpm$xy[, "est"]),

      estGammaYX_CLPM            = c(NA, lag_raw$yx[, "est"]),
      estGammaYX_RI_CLPM         = c(NA, lag_ric$yx[, "est"]),
      estGammaYX_DPM             = c(NA, lag_dpm0$yx[, "est"]),
      estGammaYX_RI_CLPM_Free    = c(NA, lag_ric_free$yx[, "est"]),
      estGammaYX_DPM_Free        = c(NA, lag_dpm_free$yx[, "est"]),
      estGammaYX_CLPM_Adj        = c(NA, lag_adj$yx[, "est"]),
      estGammaYX_CLPM_LBCA       = c(NA, lag_lbca$yx[, "est"]),
      estGammaYX_RI_CLPM_BCA     = c(NA, lag_bca_ric$yx[, "est"]),
      estGammaYX_DPM_BCA         = c(NA, lag_bca_dpm$yx[, "est"]),

      # autoregressive effects = beta
      estBetaX_CLPM            = c(NA, lag_raw$ar_x[, "est"]),
      estBetaX_RI_CLPM         = c(NA, lag_ric$ar_x[, "est"]),
      estBetaX_DPM             = c(NA, lag_dpm0$ar_x[, "est"]),
      estBetaX_RI_CLPM_Free    = c(NA, lag_ric_free$ar_x[, "est"]),
      estBetaX_DPM_Free        = c(NA, lag_dpm_free$ar_x[, "est"]),
      estBetaX_CLPM_Adj        = c(NA, lag_adj$ar_x[, "est"]),
      estBetaX_CLPM_LBCA       = c(NA, lag_lbca$ar_x[, "est"]),
      estBetaX_RI_CLPM_BCA     = c(NA, lag_bca_ric$ar_x[, "est"]),
      estBetaX_DPM_BCA         = c(NA, lag_bca_dpm$ar_x[, "est"]),

      estBetaY_CLPM            = c(NA, lag_raw$ar_y[, "est"]),
      estBetaY_RI_CLPM         = c(NA, lag_ric$ar_y[, "est"]),
      estBetaY_DPM             = c(NA, lag_dpm0$ar_y[, "est"]),
      estBetaY_RI_CLPM_Free    = c(NA, lag_ric_free$ar_y[, "est"]),
      estBetaY_DPM_Free        = c(NA, lag_dpm_free$ar_y[, "est"]),
      estBetaY_CLPM_Adj        = c(NA, lag_adj$ar_y[, "est"]),
      estBetaY_CLPM_LBCA       = c(NA, lag_lbca$ar_y[, "est"]),
      estBetaY_RI_CLPM_BCA     = c(NA, lag_bca_ric$ar_y[, "est"]),
      estBetaY_DPM_BCA         = c(NA, lag_bca_dpm$ar_y[, "est"]),

      # residual correlations from CLPM = rho
      estRho_CLPM = rho_clpm,

      # model fit failures (only non-convergence / no fit object)
      fail_CLPM           = !isTRUE(chk_clpm$converged),
      fail_RI_CLPM        = !isTRUE(chk_ric$converged),
      fail_DPM            = !isTRUE(chk_dpm0$converged),
      fail_RI_CLPM_Free   = !isTRUE(chk_ric_free$converged),
      fail_DPM_Free       = !isTRUE(chk_dpm_free$converged),
      fail_CLPM_Adj       = !isTRUE(chk_adj$converged),
      fail_CLPM_LBCA      = !isTRUE(chk_lbca$converged),
      fail_RI_CLPM_BCA    = !isTRUE(chk_bca_ric$converged),
      fail_DPM_BCA        = !isTRUE(chk_bca_dpm$converged),

      # improper solutions (converged but not proper)
      improper_CLPM           = isTRUE(chk_clpm$converged)     && !isTRUE(chk_clpm$proper),
      improper_RI_CLPM        = isTRUE(chk_ric$converged)      && !isTRUE(chk_ric$proper),
      improper_DPM            = isTRUE(chk_dpm0$converged)     && !isTRUE(chk_dpm0$proper),
      improper_RI_CLPM_Free   = isTRUE(chk_ric_free$converged) && !isTRUE(chk_ric_free$proper),
      improper_DPM_Free       = isTRUE(chk_dpm_free$converged) && !isTRUE(chk_dpm_free$proper),
      improper_CLPM_Adj       = isTRUE(chk_adj$converged)      && !isTRUE(chk_adj$proper),
      improper_CLPM_LBCA      = isTRUE(chk_lbca$converged)     && !isTRUE(chk_lbca$proper),
      improper_RI_CLPM_BCA    = isTRUE(chk_bca_ric$converged)  && !isTRUE(chk_bca_ric$proper),
      improper_DPM_BCA        = isTRUE(chk_bca_dpm$converged)  && !isTRUE(chk_bca_dpm$proper),

      # improper reasons
      improper_reason_CLPM           = rep(collapse_reasons(chk_clpm),     T),
      improper_reason_RI_CLPM        = rep(collapse_reasons(chk_ric),      T),
      improper_reason_DPM            = rep(collapse_reasons(chk_dpm0),     T),
      improper_reason_RI_CLPM_Free   = rep(collapse_reasons(chk_ric_free), T),
      improper_reason_DPM_Free       = rep(collapse_reasons(chk_dpm_free), T),
      improper_reason_CLPM_Adj       = rep(collapse_reasons(chk_adj),      T),
      improper_reason_CLPM_LBCA      = rep(collapse_reasons(chk_lbca),     T),
      improper_reason_RI_CLPM_BCA    = rep(collapse_reasons(chk_bca_ric),  T),
      improper_reason_DPM_BCA        = rep(collapse_reasons(chk_bca_dpm),  T),

      # error messages for failed models (fitting errors only)
      err_CLPM           = rep(res_clpm$err,     T),
      err_RI_CLPM        = rep(res_ric$err,      T),
      err_DPM            = rep(res_dpm0$err,     T),
      err_RI_CLPM_Free   = rep(res_ric_free$err, T),
      err_DPM_Free       = rep(res_dpm_free$err, T),
      err_CLPM_Adj       = rep(res_adj$err,      T),
      err_CLPM_LBCA      = rep(res_lbca$err,     T),
      err_RI_CLPM_BCA    = rep(res_bca_ric$err,  T),
      err_DPM_BCA        = rep(res_bca_dpm$err,  T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy[, "est"],
        lag_ric$xy[, "est"],
        lag_dpm0$xy[, "est"],
        lag_ric_free$xy[, "est"],
        lag_dpm_free$xy[, "est"],
        lag_adj$xy[, "est"],
        lag_lbca$xy[, "est"],
        lag_bca_ric$xy[, "est"],
        lag_bca_dpm$xy[, "est"]
      ))))
    )
  }

  # combine scenario outputs
  dplyr::bind_rows(out_list)
}
