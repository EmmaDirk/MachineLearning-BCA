run_one_rep_study <- function(
  rep_id,                                                 # replication index (set by outer loop)
  N,                                                      # sample size
  T,                                                      # number of waves
  k,                                                      # number of confounders
  scenarios,                                              # scenario names
  B_scenarios,                                            # named list of beta trajectories
  A,                                                      # 2×2 autoregressive + cross-lag matrix
  Psi,                                                    # k×k confounder covariance
  rho_extra,                                              # extra covariance added to X,Y each wave
  models_to_run,                                          # e.g. c("clpm","riclpm","dpm","lbca","adj")
  base_seed = 1234,                                       # base seed
  ci_level = 0.95                                         # CI level for extracted parameters
){

  # set seed for this replication
  set.seed(base_seed + rep_id)

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
    stop("Psi must be a k x k matrix with k = ", k, ".")

  # prepare output list for each scenario
  out_list <- vector("list", length(scenarios))

  # loop over the scenarios
  for (j in seq_along(scenarios)) {

    # scenario name
    scen <- scenarios[j]

    # take the beta trajectory for this scenario
    B_list <- B_scenarios[[scen]]

    # check the trajectory has length T
    if (!is.list(B_list) || length(B_list) != T)
      stop("B_scenarios[['", scen, "']] must be a list of length T.")

    # check each matrix is 2 x k
    bad_dims <- vapply(B_list, function(Bt) {
      !is.matrix(Bt) || any(dim(Bt) != c(2, k))
    }, logical(1))

    if (any(bad_dims))
      stop("One or more B matrices in scenario '", scen, "' are not 2 x k with k = ", k, ".")

    # extract mean betas
    beta_X_vec <- sapply(B_list, function(Bt) mean(Bt[1, ]))
    beta_Y_vec <- sapply(B_list, function(Bt) mean(Bt[2, ]))
    beta_vec   <- beta_X_vec

    # simulate panel data
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

    # if simulation fails, return NA rows with the correct structure
    if (is.null(df)) {

      out_list[[j]] <- data.frame(
        run      = rep(rep_id, T),
        occasion = 1:T,
        scenario = scen,

        beta     = beta_vec,
        beta_X   = beta_X_vec,
        beta_Y   = beta_Y_vec,

        estXY_CLPM      = NA,
        estXY_RI_CLPM   = NA,
        estXY_DPM       = NA,
        estXY_CLPM_Adj  = NA,
        estXY_CLPM_LBCA = NA,

        estYX_CLPM      = NA,
        estYX_RI_CLPM   = NA,
        estYX_DPM       = NA,
        estYX_CLPM_Adj  = NA,
        estYX_CLPM_LBCA = NA,

        estA_CLPM      = NA,
        estA_RI_CLPM   = NA,
        estA_DPM       = NA,
        estA_CLPM_Adj  = NA,
        estA_CLPM_LBCA = NA,

        estAY_CLPM      = NA,
        estAY_RI_CLPM   = NA,
        estAY_DPM       = NA,
        estAY_CLPM_Adj  = NA,
        estAY_CLPM_LBCA = NA,

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
    res_clpm <- if ("clpm"   %in% models_to_run) safe_fit_clpm(model_clpm, df)           else list(fit = NULL, err = NA_character_)
    res_ric  <- if ("riclpm" %in% models_to_run) safe_fit_riclpm(model_riclpm, df)       else list(fit = NULL, err = NA_character_)
    res_dpm0 <- if ("dpm"    %in% models_to_run) safe_fit_dpm(model_dpm, df)             else list(fit = NULL, err = NA_character_)
    res_adj  <- if ("adj"    %in% models_to_run) safe_fit_clpm_C(model_clpm_with_Cs, df) else list(fit = NULL, err = NA_character_)
    res_lbca <- if ("lbca"   %in% models_to_run) safe_fit_clpm_resid(model_clpm, df)     else list(fit = NULL, err = NA_character_)

    # pull fitted objects
    fit_clpm_raw <- res_clpm$fit
    fit_ric      <- res_ric$fit
    fit_dpm0     <- res_dpm0$fit
    fit_adj      <- res_adj$fit
    fit_lbca     <- res_lbca$fit

    # extract lagged parameters
    lag_raw  <- extract_lagged_parameters(fit_clpm_raw, T, "clpm",   ci_level = ci_level)
    lag_ric  <- extract_lagged_parameters(fit_ric,       T, "riclpm", ci_level = ci_level)
    lag_dpm0 <- extract_lagged_parameters(fit_dpm0,      T, "dpm",    ci_level = ci_level)
    lag_adj  <- extract_lagged_parameters(fit_adj,       T, "clpm",   ci_level = ci_level)
    lag_lbca <- extract_lagged_parameters(fit_lbca,      T, "clpm",   ci_level = ci_level)

    # extract residual correlations from CLPM
    rho_clpm <- extract_rho_vec(fit_clpm_raw, T, "clpm")

    # assemble output rows
    out_list[[j]] <- data.frame(

      run      = rep(rep_id, T),
      occasion = 1:T,
      scenario = scen,

      beta     = beta_vec,
      beta_X   = beta_X_vec,
      beta_Y   = beta_Y_vec,

      estXY_CLPM      = c(NA, lag_raw$xy[, "est"]),
      estXY_RI_CLPM   = c(NA, lag_ric$xy[, "est"]),
      estXY_DPM       = c(NA, lag_dpm0$xy[, "est"]),
      estXY_CLPM_Adj  = c(NA, lag_adj$xy[, "est"]),
      estXY_CLPM_LBCA = c(NA, lag_lbca$xy[, "est"]),

      estYX_CLPM      = c(NA, lag_raw$yx[, "est"]),
      estYX_RI_CLPM   = c(NA, lag_ric$yx[, "est"]),
      estYX_DPM       = c(NA, lag_dpm0$yx[, "est"]),
      estYX_CLPM_Adj  = c(NA, lag_adj$yx[, "est"]),
      estYX_CLPM_LBCA = c(NA, lag_lbca$yx[, "est"]),

      estA_CLPM      = c(NA, lag_raw$ar_x[, "est"]),
      estA_RI_CLPM   = c(NA, lag_ric$ar_x[, "est"]),
      estA_DPM       = c(NA, lag_dpm0$ar_x[, "est"]),
      estA_CLPM_Adj  = c(NA, lag_adj$ar_x[, "est"]),
      estA_CLPM_LBCA = c(NA, lag_lbca$ar_x[, "est"]),

      estAY_CLPM      = c(NA, lag_raw$ar_y[, "est"]),
      estAY_RI_CLPM   = c(NA, lag_ric$ar_y[, "est"]),
      estAY_DPM       = c(NA, lag_dpm0$ar_y[, "est"]),
      estAY_CLPM_Adj  = c(NA, lag_adj$ar_y[, "est"]),
      estAY_CLPM_LBCA = c(NA, lag_lbca$ar_y[, "est"]),

      estRho_CLPM = rho_clpm,

      fail_CLPM      = is.null(fit_clpm_raw),
      fail_RI_CLPM   = is.null(fit_ric),
      fail_DPM       = is.null(fit_dpm0),
      fail_CLPM_Adj  = is.null(fit_adj),
      fail_CLPM_LBCA = is.null(fit_lbca),

      err_CLPM      = rep(res_clpm$err,  T),
      err_RI_CLPM   = rep(res_ric$err,   T),
      err_DPM       = rep(res_dpm0$err,  T),
      err_CLPM_Adj  = rep(res_adj$err,   T),
      err_CLPM_LBCA = rep(res_lbca$err,  T),

      is_na_run = as.integer(all(is.na(c(
        lag_raw$xy[, "est"],
        lag_ric$xy[, "est"],
        lag_dpm0$xy[, "est"],
        lag_adj$xy[, "est"],
        lag_lbca$xy[, "est"]
      ))))
    )
  }

  dplyr::bind_rows(out_list)
}
