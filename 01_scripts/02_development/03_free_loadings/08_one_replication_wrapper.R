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
# 6) It calls the appropriate function to extract parameter estimates and confidence intervals.
# 7) It stores the results in a list.
# ------------------------------------------------------------------------------------------------------------

run_one_rep_study <- function(
  rep_id,                                                 # replication index (set by outer loop)
  N,                                                      # sample size
  T,                                                      # number of waves
  k,                                                      # number of confounders
  scenarios,                                              # scenario names
  D_scenarios,                                            # named list of delta trajectories
  A,                                                      # 2×2 autoregressive (beta) + cross-lagged (gamma) matrix
  Psi,                                                    # k×k confounder covariance
  rho_extra = 0.15,                                       # extra covariance added to X,Y at time t
  models_to_run = c("clpm","riclpm","riclpm_free_RI_loadings","dpm","dpm_free_loadings","lbca","adj"),  # default is all models
  base_seed = 1234,                                       # base seed
  ci_level = 0.95                                         # CI level for extracted parameters
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

        # delta (confounder effects) 
        delta    = delta_vec,
        delta_X  = delta_X_vec,
        delta_Y  = delta_Y_vec,

        # cross-lagged effects = gamma
        estGammaXY_CLPM                    = NA,
        estGammaXY_RI_CLPM                 = NA,
        estGammaXY_RI_CLPM_free_RI_loadings = NA,
        estGammaXY_DPM                     = NA,
        estGammaXY_DPM_free_loadings       = NA,
        estGammaXY_CLPM_Adj                = NA,
        estGammaXY_CLPM_LBCA               = NA,

        estGammaYX_CLPM                    = NA,
        estGammaYX_RI_CLPM                 = NA,
        estGammaYX_RI_CLPM_free_RI_loadings = NA,
        estGammaYX_DPM                     = NA,
        estGammaYX_DPM_free_loadings       = NA,
        estGammaYX_CLPM_Adj                = NA,
        estGammaYX_CLPM_LBCA               = NA,

        # autoregressive effects = beta
        estBetaX_CLPM                    = NA,
        estBetaX_RI_CLPM                 = NA,
        estBetaX_RI_CLPM_free_RI_loadings = NA,
        estBetaX_DPM                     = NA,
        estBetaX_DPM_free_loadings       = NA,
        estBetaX_CLPM_Adj                = NA,
        estBetaX_CLPM_LBCA               = NA,

        estBetaY_CLPM                    = NA,
        estBetaY_RI_CLPM                 = NA,
        estBetaY_RI_CLPM_free_RI_loadings = NA,
        estBetaY_DPM                     = NA,
        estBetaY_DPM_free_loadings       = NA,
        estBetaY_CLPM_Adj                = NA,
        estBetaY_CLPM_LBCA               = NA,

        # residual correlations from CLPM = rho
        estRho_CLPM = NA,
        
        # model fit failures
        fail_CLPM                    = TRUE,
        fail_RI_CLPM                 = TRUE,
        fail_RI_CLPM_free_RI_loadings = TRUE,
        fail_DPM                     = TRUE,
        fail_DPM_free_loadings       = TRUE,
        fail_CLPM_Adj                = TRUE,
        fail_CLPM_LBCA               = TRUE,

        # error messages if sim failed
        err_CLPM                    = "sim failed",
        err_RI_CLPM                 = "sim failed",
        err_RI_CLPM_free_RI_loadings = "sim failed",
        err_DPM                     = "sim failed",
        err_DPM_free_loadings       = "sim failed",
        err_CLPM_Adj                = "sim failed",
        err_CLPM_LBCA               = "sim failed",
        
        # NA run marker
        is_na_run = 1L
      )

      next
    }

    # build model strings
    model_clpm                 <- build_clpm(T)
    model_riclpm               <- build_riclpm(T)
    model_riclpm_free_RI       <- build_riclpm_free_RI_loadings(T)
    model_dpm                  <- build_dpm(T)
    model_dpm_free_loadings    <- build_dpm_free_loadings(T)
    model_clpm_with_Cs         <- build_clpm_with_Cs(T, k)

    # fit the models safely
    res_clpm     <- if ("clpm"                     %in% models_to_run) safe_fit_clpm(model_clpm, df)                         else list(fit = NULL, err = NA_character_)
    res_ric      <- if ("riclpm"                   %in% models_to_run) safe_fit_riclpm(model_riclpm, df)                     else list(fit = NULL, err = NA_character_)
    res_ric_free <- if ("riclpm_free_RI_loadings"  %in% models_to_run) safe_fit_riclpm_free_RI_loadings(model_riclpm_free_RI, df) else list(fit = NULL, err = NA_character_)
    res_dpm0     <- if ("dpm"                      %in% models_to_run) safe_fit_dpm(model_dpm, df)                           else list(fit = NULL, err = NA_character_)
    res_dpm_free <- if ("dpm_free_loadings"        %in% models_to_run) safe_fit_dpm_free_loadings(model_dpm_free_loadings, df) else list(fit = NULL, err = NA_character_)
    res_adj      <- if ("adj"                      %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df)               else list(fit = NULL, err = NA_character_)
    res_lbca     <- if ("lbca"                     %in% models_to_run) safe_fit_clpm_resid(model_clpm, df)                   else list(fit = NULL, err = NA_character_)

    # pull fitted objects
    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_ric_free <- res_ric_free$fit
    fit_dpm0     <- res_dpm0$fit
    fit_dpm_free <- res_dpm_free$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit

    # extract lagged parameters
    lag_raw      <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",                    ci_level = ci_level)
    lag_ric      <- extract_lagged_parameters(fit_ric,       T, "riclpm",                  ci_level = ci_level)
    lag_ric_free <- extract_lagged_parameters(fit_ric_free,  T, "riclpm_free_RI_loadings", ci_level = ci_level)
    lag_dpm0     <- extract_lagged_parameters(fit_dpm0,      T, "dpm",                     ci_level = ci_level)
    lag_dpm_free <- extract_lagged_parameters(fit_dpm_free,  T, "dpm_free_loadings",       ci_level = ci_level)
    lag_adj      <- extract_lagged_parameters(fit_adj,       T, "clpm",                    ci_level = ci_level)
    lag_lbca     <- extract_lagged_parameters(fit_lbca,      T, "clpm",                    ci_level = ci_level)

    # extract residual correlations from CLPM
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm")

    # assemble output rows
    out_list[[j]] <- data.frame(

      # run info
      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      # delta (confounder effects) 
      delta    = delta_vec,
      delta_X  = delta_X_vec,
      delta_Y  = delta_Y_vec,

      # cross-lagged effects = gamma
      estGammaXY_CLPM                    = c(NA, lag_raw$xy[, "est"]),
      estGammaXY_RI_CLPM                 = c(NA, lag_ric$xy[, "est"]),
      estGammaXY_RI_CLPM_free_RI_loadings = c(NA, lag_ric_free$xy[, "est"]),
      estGammaXY_DPM                     = c(NA, lag_dpm0$xy[, "est"]),
      estGammaXY_DPM_free_loadings       = c(NA, lag_dpm_free$xy[, "est"]),
      estGammaXY_CLPM_Adj                = c(NA, lag_adj$xy[, "est"]),
      estGammaXY_CLPM_LBCA               = c(NA, lag_lbca$xy[, "est"]),

      estGammaYX_CLPM                    = c(NA, lag_raw$yx[, "est"]),
      estGammaYX_RI_CLPM                 = c(NA, lag_ric$yx[, "est"]),
      estGammaYX_RI_CLPM_free_RI_loadings = c(NA, lag_ric_free$yx[, "est"]),
      estGammaYX_DPM                     = c(NA, lag_dpm0$yx[, "est"]),
      estGammaYX_DPM_free_loadings       = c(NA, lag_dpm_free$yx[, "est"]),
      estGammaYX_CLPM_Adj                = c(NA, lag_adj$yx[, "est"]),
      estGammaYX_CLPM_LBCA               = c(NA, lag_lbca$yx[, "est"]),

      # autoregressive effects = beta
      estBetaX_CLPM                    = c(NA, lag_raw$ar_x[, "est"]),
      estBetaX_RI_CLPM                 = c(NA, lag_ric$ar_x[, "est"]),
      estBetaX_RI_CLPM_free_RI_loadings = c(NA, lag_ric_free$ar_x[, "est"]),
      estBetaX_DPM                     = c(NA, lag_dpm0$ar_x[, "est"]),
      estBetaX_DPM_free_loadings       = c(NA, lag_dpm_free$ar_x[, "est"]),
      estBetaX_CLPM_Adj                = c(NA, lag_adj$ar_x[, "est"]),
      estBetaX_CLPM_LBCA               = c(NA, lag_lbca$ar_x[, "est"]),

      estBetaY_CLPM                    = c(NA, lag_raw$ar_y[, "est"]),
      estBetaY_RI_CLPM                 = c(NA, lag_ric$ar_y[, "est"]),
      estBetaY_RI_CLPM_free_RI_loadings = c(NA, lag_ric_free$ar_y[, "est"]),
      estBetaY_DPM                     = c(NA, lag_dpm0$ar_y[, "est"]),
      estBetaY_DPM_free_loadings       = c(NA, lag_dpm_free$ar_y[, "est"]),
      estBetaY_CLPM_Adj                = c(NA, lag_adj$ar_y[, "est"]),
      estBetaY_CLPM_LBCA               = c(NA, lag_lbca$ar_y[, "est"]),

      # residual correlations from CLPM = rho
      estRho_CLPM = rho_clpm,

      # model fit failures
      fail_CLPM                    = is.null(fit_clpm_raw),
      fail_RI_CLPM                 = is.null(fit_ric),
      fail_RI_CLPM_free_RI_loadings = is.null(fit_ric_free),
      fail_DPM                     = is.null(fit_dpm0),
      fail_DPM_free_loadings       = is.null(fit_dpm_free),
      fail_CLPM_Adj                = is.null(fit_adj),
      fail_CLPM_LBCA               = is.null(fit_lbca),

      # error messages for failed models
      err_CLPM                    = rep(res_clpm$err,     T),
      err_RI_CLPM                 = rep(res_ric$err,      T),
      err_RI_CLPM_free_RI_loadings = rep(res_ric_free$err, T),
      err_DPM                     = rep(res_dpm0$err,     T),
      err_DPM_free_loadings       = rep(res_dpm_free$err, T),
      err_CLPM_Adj                = rep(res_adj$err,      T),
      err_CLPM_LBCA               = rep(res_lbca$err,     T),

      # NA run marker
      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy[, "est"],
        lag_ric$xy[, "est"],
        lag_ric_free$xy[, "est"],
        lag_dpm0$xy[, "est"],
        lag_dpm_free$xy[, "est"],
        lag_adj$xy[, "est"],
        lag_lbca$xy[, "est"]
      ))))
    )
  }

  # combine scenario outputs
  dplyr::bind_rows(out_list)
}
