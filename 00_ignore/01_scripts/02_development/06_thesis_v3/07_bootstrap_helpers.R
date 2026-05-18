# 07_bootstrap_helpers.R
# These helpers handle the parts of the workflow that are not the main SEM fit itself:
# - classifying runs as success / non-converged / improper
# - bootstrap-based standard errors for two-stage procedures
#
# The bootstrap logic is now multi-model aware:
# - each bootstrap draw samples the data only once
# - each unique stage-1 recipe is prepared only once inside that draw
# - every compatible SEM is then fit on that same prepared bootstrap sample
#
# Flag coding:
# - 0 = successful and proper run
# - 1 = failed or non-converged run
# - 2 = converged but improper solution
# -------------------------------------------------------------------------------------------------

# classify a fitted lavaan object into the requested flag coding
classify_fit_flag <- function(fit) {

  # completely failed fit counts as non-convergence
  if (is.null(fit)) {
    return(1L)
  }

  # check convergence first
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) {
    return(1L)
  }

  # lavaan already provides a broad admissibility check
  post_ok <- tryCatch(lavaan::lavInspect(fit, "post.check"), error = function(e) TRUE)
  if (!isTRUE(post_ok)) {
    return(2L)
  }

  # extra guard: catch negative variances if they appear in the parameter table
  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)
  if (!is.null(pe)) {
    var_rows <- pe$op == "~~" & pe$lhs == pe$rhs
    if (any(pe$est[var_rows] < 0, na.rm = TRUE)) {
      return(2L)
    }
  }

  # otherwise it is a successful run
  0L
}


# convert one integer flag into three proportions that sum to 1
flag_to_props <- function(flag) {

  out <- c(flag0 = 0, flag1 = 0, flag2 = 0)

  if (is.na(flag) || !(flag %in% c(0L, 1L, 2L, 0, 1, 2))) {
    return(as.list(out))
  }

  out[paste0("flag", as.integer(flag))] <- 1

  as.list(out)
}


# helper to turn a variance name into a readable diagnosis
describe_negative_variance_target <- function(name) {

  if (is.na(name) || !nzchar(name)) {
    return("negative variance")
  }

  # RI-CLPM random intercepts
  if (grepl("^ri[xy]$", name)) {
    return(sprintf("negative random intercept variance (%s)", name))
  }

  # DPM accumulating factors
  if (grepl("^F[XY]$", name)) {
    return(sprintf("negative accumulating factor variance (%s)", name))
  }

  # RI-CLPM within-person factors
  if (grepl("^w[xy][0-9]+$", name)) {
    return(sprintf("negative within-person latent variance (%s)", name))
  }

  # observed variables / residual variances
  if (grepl("^[xy][0-9]+$", name)) {
    return(sprintf("negative observed residual variance (%s)", name))
  }

  # generic latent or observed variance
  sprintf("negative variance (%s)", name)
}


# helper to detect whether a covariance matrix is not positive semidefinite
matrix_not_psd <- function(S, tol = 1e-10) {

  if (is.null(S)) {
    return(FALSE)
  }

  if (!is.matrix(S)) {
    S <- tryCatch(as.matrix(S), error = function(e) NULL)
  }

  if (is.null(S)) {
    return(FALSE)
  }

  if (nrow(S) == 0 || ncol(S) == 0) {
    return(FALSE)
  }

  ev <- tryCatch(
    eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  if (all(is.na(ev))) {
    return(FALSE)
  }

  any(ev < -tol, na.rm = TRUE)
}


# diagnose the main reason why a converged fit is improper
diagnose_improper_fit <- function(fit, tol = 1e-10) {

  # failed fit
  if (is.null(fit)) {
    return(NA_character_)
  }

  # non-converged fit
  converged <- tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) {
    return(NA_character_)
  }

  # if lavaan thinks the fit is admissible, do not attach a reason
  post_ok <- tryCatch(lavaan::lavInspect(fit, "post.check"), error = function(e) TRUE)

  pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)

  # 1) first and most interpretable case: explicit negative variances
  if (!is.null(pe)) {

    var_rows <- pe$op == "~~" & pe$lhs == pe$rhs & !is.na(pe$est)

    if (any(var_rows)) {

      neg_var_rows <- which(var_rows & pe$est < -tol)

      if (length(neg_var_rows) > 0) {

        # choose the most negative variance as the main culprit
        i <- neg_var_rows[which.min(pe$est[neg_var_rows])]
        return(describe_negative_variance_target(pe$lhs[i]))
      }
    }
  }

  # 2) check the latent covariance matrix
  cov_lv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
  if (matrix_not_psd(cov_lv, tol = tol)) {
    return("latent covariance matrix not positive semidefinite")
  }

  # 3) check the observed residual covariance matrix
  theta <- tryCatch(lavaan::lavInspect(fit, "theta"), error = function(e) NULL)
  if (matrix_not_psd(theta, tol = tol)) {
    return("observed residual covariance matrix not positive semidefinite")
  }

  # 4) if lavaan flagged the solution as improper but we did not localize it more precisely
  if (!isTRUE(post_ok)) {
    return("post.check failed: unspecified improper solution")
  }

  # otherwise there is no improper-fit reason to report
  NA_character_
}


# convert fit status into a compact bootstrap issue label
diagnose_bootstrap_issue <- function(fit, tol = 1e-10) {

  flag <- classify_fit_flag(fit)

  if (flag == 0L) {
    return("proper")
  }

  if (flag == 1L) {
    return("nonconverged_or_failed")
  }

  reason <- diagnose_improper_fit(fit, tol = tol)

  if (is.na(reason) || !nzchar(reason)) {
    return("improper_unspecified")
  }

  reason
}


# determine whether the chosen pipeline needs bootstrap-based standard errors
uses_bootstrap_se <- function(residualizer) {

  residualizer %in% c("linear", "xgb")
}


# make one empty bootstrap summary for one model
make_empty_bootstrap_summary <- function(T) {
  list(
    ARX = rep(NA_real_, T),
    ARY = rep(NA_real_, T),
    CXY = rep(NA_real_, T),
    CYX = rep(NA_real_, T),
    mean_boot_mse_x = rep(NA_real_, T),
    sd_boot_mse_x = rep(NA_real_, T),
    mean_boot_r2_x = rep(NA_real_, T),
    sd_boot_r2_x = rep(NA_real_, T),
    mean_boot_mse_y = rep(NA_real_, T),
    sd_boot_mse_y = rep(NA_real_, T),
    mean_boot_r2_y = rep(NA_real_, T),
    sd_boot_r2_y = rep(NA_real_, T),
    bootstrap_prop_success = NA_real_,
    flag0 = NA_real_,
    flag1 = NA_real_,
    flag2 = NA_real_,
    bootstrap_issue_vector = NA_character_
  )
}


# summarize one bootstrap metric matrix column-wise while keeping all-NA columns as NA
summarise_boot_metric <- function(M, fun = mean) {

  apply(M, 2, function(z) {
    if (all(is.na(z))) {
      return(NA_real_)
    }
    as.numeric(fun(z, na.rm = TRUE))
  })
}


# bootstrap a whole set of models while sharing samples and stage-1 preparation
bootstrap_model_set <- function(
    df,
    T,
    k,
    model_specs,
    stage1_groups,
    seed = NULL
) {

  # keep only the model specs that actually require bootstrap-based SEs
  bootstrap_specs <- Filter(function(x) uses_bootstrap_se(x$residualizer) && x$bootstrap_B >= 2L,
                            model_specs)

  # if no model needs bootstrap, return an empty named list
  if (length(bootstrap_specs) == 0) {
    return(setNames(vector("list", 0), character(0)))
  }

  # if no bootstrap seed was supplied, generate one automatically
  if (is.null(seed)) {
    max_seed <- max(1L, .Machine$integer.max - 100000L)
    seed <- as.integer(sample.int(max_seed, size = 1))
  } else {
    seed <- as.integer(seed[1])
  }

  # set the bootstrap seed once
  set.seed(seed)

  # ensure that each original row carries a stable id through the bootstrap
  if (!(".id_orig" %in% names(df))) {
    df$.id_orig <- seq_len(nrow(df))
  }

  # maximum number of draws required by any bootstrap model
  max_B <- max(vapply(bootstrap_specs, function(x) x$bootstrap_B, integer(1)))

  # storage for every model separately
  store <- setNames(vector("list", length(bootstrap_specs)),
                    vapply(bootstrap_specs, function(x) x$name, character(1)))

  for (spec in bootstrap_specs) {

    B_spec <- spec$bootstrap_B

    store[[spec$name]] <- list(
      B = B_spec,
      ARX = matrix(NA_real_, nrow = B_spec, ncol = T),
      ARY = matrix(NA_real_, nrow = B_spec, ncol = T),
      CXY = matrix(NA_real_, nrow = B_spec, ncol = T),
      CYX = matrix(NA_real_, nrow = B_spec, ncol = T),
      mse_x = matrix(NA_real_, nrow = B_spec, ncol = T),
      r2_x = matrix(NA_real_, nrow = B_spec, ncol = T),
      mse_y = matrix(NA_real_, nrow = B_spec, ncol = T),
      r2_y = matrix(NA_real_, nrow = B_spec, ncol = T),
      boot_flag = rep(NA_integer_, B_spec),
      bootstrap_issue_vector = rep(NA_character_, B_spec)
    )
  }

  # only stage-1 groups that are needed by bootstrap models matter here
  bootstrap_group_ids <- sort(unique(vapply(bootstrap_specs, function(x) x$stage1_group_id, integer(1))))
  bootstrap_groups <- Filter(function(g) g$stage1_group_id %in% bootstrap_group_ids, stage1_groups)

  # bootstrap the shared pipeline
  for (b in seq_len(max_B)) {

    # sample rows with replacement once
    idx <- sample.int(n = nrow(df), size = nrow(df), replace = TRUE)
    df_b <- df[idx, , drop = FALSE]

    # prepare every required stage-1 group exactly once inside this draw
    prepared_by_group <- list()

    for (g in seq_along(bootstrap_groups)) {

      group_obj <- bootstrap_groups[[g]]
      proto <- group_obj$prototype

      residualizer_args_b <- proto$residualizer_args

      if (is.null(residualizer_args_b)) {
        residualizer_args_b <- list()
      }

      # vary the stage-1 seed across draws, but keep it shared within the group
      residualizer_args_b$seed <- as.integer(seed + 1000L * g + b)

      prepared_by_group[[as.character(group_obj$stage1_group_id)]] <- prepare_analysis_data(
        df = df_b,
        k = k,
        residualizer = proto$residualizer,
        confounder_order = proto$confounder_order,
        exclude = proto$exclude,
        xgb_tuning = proto$xgb_tuning,
        residualizer_args = residualizer_args_b
      )
    }

    # fit every bootstrap-using SEM on its already prepared data
    for (spec in bootstrap_specs) {

      if (b > spec$bootstrap_B) {
        next
      }

      prep <- prepared_by_group[[as.character(spec$stage1_group_id)]]

      # store the stage-1 ML diagnostics for this bootstrap draw
      if (!is.null(prep$ml_metrics)) {
        boot_ml <- standardize_ml_metric_frame(prep$ml_metrics, T = T)
        store[[spec$name]]$mse_x[b, ] <- boot_ml$mse_x
        store[[spec$name]]$r2_x[b, ] <- boot_ml$r2_x
        store[[spec$name]]$mse_y[b, ] <- boot_ml$mse_y
        store[[spec$name]]$r2_y[b, ] <- boot_ml$r2_y
      }

      if (is.null(prep$data)) {
        fit_b <- list(fit = NULL)
      } else {
        fit_b <- fit_sem_on_prepared_data(
          df_prepared = prep$data,
          T = T,
          k = k,
          residualizer = spec$residualizer,
          sem_model = spec$sem_model,
          confounder_order = spec$confounder_order,
          exclude = spec$exclude,
          free_loadings = spec$free_loadings
        )
      }

      # classify every bootstrap fit, including outright failures
      store[[spec$name]]$boot_flag[b] <- classify_fit_flag(fit_b$fit)

      # store a readable bootstrap issue label
      store[[spec$name]]$bootstrap_issue_vector[b] <- diagnose_bootstrap_issue(fit_b$fit)

      # skip failed bootstrap fits for estimate extraction
      if (is.null(fit_b$fit)) {
        next
      }

      # extract lagged estimates from the bootstrap fit
      lag_b <- extract_lagged_estimates(
        fit = fit_b$fit,
        T = T,
        model_type = spec$sem_model
      )

      store[[spec$name]]$ARX[b, ] <- lag_b$ARX
      store[[spec$name]]$ARY[b, ] <- lag_b$ARY
      store[[spec$name]]$CXY[b, ] <- lag_b$CXY
      store[[spec$name]]$CYX[b, ] <- lag_b$CYX
    }
  }

  # collapse bootstrap storage into the final summaries
  out <- setNames(vector("list", length(bootstrap_specs)), names(store))

  for (nm in names(store)) {

    obj <- store[[nm]]

    out[[nm]] <- list(
      ARX = summarise_boot_metric(obj$ARX, fun = stats::sd),
      ARY = summarise_boot_metric(obj$ARY, fun = stats::sd),
      CXY = summarise_boot_metric(obj$CXY, fun = stats::sd),
      CYX = summarise_boot_metric(obj$CYX, fun = stats::sd),
      mean_boot_mse_x = summarise_boot_metric(obj$mse_x, fun = mean),
      sd_boot_mse_x = summarise_boot_metric(obj$mse_x, fun = stats::sd),
      mean_boot_r2_x = summarise_boot_metric(obj$r2_x, fun = mean),
      sd_boot_r2_x = summarise_boot_metric(obj$r2_x, fun = stats::sd),
      mean_boot_mse_y = summarise_boot_metric(obj$mse_y, fun = mean),
      sd_boot_mse_y = summarise_boot_metric(obj$mse_y, fun = stats::sd),
      mean_boot_r2_y = summarise_boot_metric(obj$r2_y, fun = mean),
      sd_boot_r2_y = summarise_boot_metric(obj$r2_y, fun = stats::sd),
      bootstrap_prop_success = mean(obj$boot_flag == 0L, na.rm = TRUE),
      flag0 = mean(obj$boot_flag == 0L, na.rm = TRUE),
      flag1 = mean(obj$boot_flag == 1L, na.rm = TRUE),
      flag2 = mean(obj$boot_flag == 2L, na.rm = TRUE),
      bootstrap_issue_vector = obj$bootstrap_issue_vector
    )
  }

  out
}


# convenience wrapper for bootstrapping a single model through the shared multi-model engine
bootstrap_pipeline_se <- function(
    df,
    T,
    k,
    residualizer,
    sem_model,
    confounder_order,
    exclude,
    free_loadings,
    xgb_tuning,
    residualizer_args,
    B,
    seed = NULL
) {

  # if B < 2, return the classic empty structure
  if (is.null(B) || B < 2) {
    return(make_empty_bootstrap_summary(T))
  }

  spec <- make_model_spec(
    name = "bootstrap_target",
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    bootstrap_B = B,
    xgb_tuning = xgb_tuning,
    tune_xgb = FALSE,
    xgb_tune_args = list(),
    residualizer_args = residualizer_args
  )

  specs <- normalize_model_spec_list(list(spec))
  specs <- assign_xgb_tuning_group_ids(specs)
  specs <- assign_stage1_group_ids(specs)
  stage1_groups <- build_stage1_groups(specs)

  out <- bootstrap_model_set(
    df = df,
    T = T,
    k = k,
    model_specs = specs,
    stage1_groups = stage1_groups,
    seed = seed
  )

  out[[spec$name]]
}