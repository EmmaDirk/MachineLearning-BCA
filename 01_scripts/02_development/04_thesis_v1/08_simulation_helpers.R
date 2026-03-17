# This script contains a few small helpers for simulation studies.
#
# The philosophy is the same as in the rest of the engine:
# - keep the helpers small
# - keep them easy to test
# - return tidy objects that are easy to bind together later
#
# The important idea is that the modelling functions are now the real workhorses.
# The simulation wrapper simply calls them repeatedly.
# That means the simulation code does not need a giant model-specific pipeline.
# ------------------------------------------------------------------------------------------

# helper: fit several models to one data set
# purpose:
# model_specs is a named list.
# each entry must contain:
# - fun  : a fitting function, e.g. fit_clpm
# - args : a list of arguments passed to that function
#
# example:
# model_specs <- list(
#   clpm = list(
#     fun = fit_clpm,
#     args = list(T = 5, residualiser = "none")
#   ),
#   riclpm = list(
#     fun = fit_riclpm,
#     args = list(T = 5, residualiser = "linear", bootstrap_R = 199)
#   )
# )
fit_model_set <- function(df, model_specs) {

  if (!is.list(model_specs) || length(model_specs) == 0) {
    stop("'model_specs' must be a non-empty list.")
  }

  out <- vector("list", length(model_specs))
  names(out) <- names(model_specs)

  for (nm in names(model_specs)) {

    spec <- model_specs[[nm]]

    if (!is.list(spec) || is.null(spec$fun) || is.null(spec$args)) {
      stop("Each entry of 'model_specs' must be a list with elements 'fun' and 'args'.")
    }

    fit_obj <- tryCatch(
      do.call(spec$fun, c(list(df = df), spec$args)),
      error = function(e) e
    )

    # if the top-level call failed before a normal fit object could be created,
    # we still return a very small placeholder.
    if (inherits(fit_obj, "error")) {
      out[[nm]] <- list(
        model_type = nm,
        residualiser = NA_character_,
        converged = FALSE,
        proper = FALSE,
        bic = NA_real_,
        se_type = NA_character_,
        parameters = data.frame(),
        preprocess_info = NULL,
        error = conditionMessage(fit_obj)
      )
    } else {
      out[[nm]] <- fit_obj
    }
  }

  out
}

# helper: convert one fit result to a tidy data frame
# purpose:
# this is the object that we can row-bind over models, scenarios, and replications.
format_fit_for_simulation <- function(
  fit_object,
  rep_id,
  scenario,
  model_name
) {

  # if the fit did not even produce a parameter table, return one empty row
  if (is.null(fit_object$parameters) || nrow(fit_object$parameters) == 0) {
    return(data.frame(
      rep_id = rep_id,
      scenario = scenario,
      model_name = model_name,
      model_type = fit_object$model_type %||% NA_character_,
      residualiser = fit_object$residualiser %||% NA_character_,
      converged = fit_object$converged %||% FALSE,
      proper = fit_object$proper %||% FALSE,
      bic = fit_object$bic %||% NA_real_,
      se_type = fit_object$se_type %||% NA_character_,
      param_id = NA_character_,
      lhs = NA_character_,
      rhs = NA_character_,
      path_type = NA_character_,
      wave_from = NA_integer_,
      wave_to = NA_integer_,
      est = NA_real_,
      se = NA_real_,
      se_model = NA_real_,
      se_boot = NA_real_,
      error = fit_object$error %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  params <- fit_object$parameters

  data.frame(
    rep_id = rep_id,
    scenario = scenario,
    model_name = model_name,
    model_type = fit_object$model_type,
    residualiser = fit_object$residualiser,
    converged = fit_object$converged,
    proper = fit_object$proper,
    bic = fit_object$bic,
    se_type = fit_object$se_type,
    param_id = params$param_id,
    lhs = params$lhs,
    rhs = params$rhs,
    path_type = params$path_type,
    wave_from = params$wave_from,
    wave_to = params$wave_to,
    est = params$est,
    se = params$se,
    se_model = params$se_model,
    se_boot = params$se_boot,
    error = fit_object$error,
    stringsAsFactors = FALSE
  )
}

# helper: add the true lagged effects to a tidy result table
# purpose:
# in your simulation study, this makes bias and RMSE calculations easier later.
add_true_lagged_effects <- function(results_df, Phi) {

  if (!is.matrix(Phi) || any(dim(Phi) != c(2, 2))) {
    stop("'Phi' must be a 2 x 2 matrix.")
  }

  results_df$true_value <- NA_real_

  results_df$true_value[results_df$path_type == "ar_x"] <- Phi[1, 1]
  results_df$true_value[results_df$path_type == "cl_y_to_x"] <- Phi[1, 2]
  results_df$true_value[results_df$path_type == "cl_x_to_y"] <- Phi[2, 1]
  results_df$true_value[results_df$path_type == "ar_y"] <- Phi[2, 2]

  results_df
}

# helper: run one replication over one or more delta scenarios
# purpose:
# for each scenario we:
# 1) simulate one data set
# 2) fit all requested models
# 3) stack the tidy outputs
run_one_replication <- function(
  rep_id,
  N,
  T,
  Delta_scenarios,
  Phi,
  Omega11,
  Sigma,
  model_specs,
  base_seed = 1234
) {

  set.seed(base_seed + rep_id)

  if (!is.list(Delta_scenarios) || length(Delta_scenarios) == 0) {
    stop("'Delta_scenarios' must be a non-empty named list.")
  }

  if (is.null(names(Delta_scenarios)) || any(names(Delta_scenarios) == "")) {
    stop("'Delta_scenarios' must be a named list.")
  }

  out_list <- vector("list", length(Delta_scenarios))

  scen_names <- names(Delta_scenarios)

  for (j in seq_along(Delta_scenarios)) {

    scen <- scen_names[j]
    Delta_list <- Delta_scenarios[[j]]

    # simulate one data set for this scenario
    df <- tryCatch(
      simulate_panel_data(
        N = N,
        T = T,
        Phi = Phi,
        Delta_list = Delta_list,
        Omega11 = Omega11,
        Sigma = Sigma,
        return_full = FALSE,
        seed = base_seed + rep_id + j
      ),
      error = function(e) e
    )

    # if the simulation itself failed, return one small row
    if (inherits(df, "error")) {
      out_list[[j]] <- data.frame(
        rep_id = rep_id,
        scenario = scen,
        model_name = NA_character_,
        model_type = NA_character_,
        residualiser = NA_character_,
        converged = FALSE,
        proper = FALSE,
        bic = NA_real_,
        se_type = NA_character_,
        param_id = NA_character_,
        lhs = NA_character_,
        rhs = NA_character_,
        path_type = NA_character_,
        wave_from = NA_integer_,
        wave_to = NA_integer_,
        est = NA_real_,
        se = NA_real_,
        se_model = NA_real_,
        se_boot = NA_real_,
        error = conditionMessage(df),
        stringsAsFactors = FALSE
      )
      next
    }

    fits <- fit_model_set(df = df, model_specs = model_specs)

    scen_rows <- do.call(rbind, lapply(names(fits), function(model_name) {
      format_fit_for_simulation(
        fit_object = fits[[model_name]],
        rep_id = rep_id,
        scenario = scen,
        model_name = model_name
      )
    }))

    out_list[[j]] <- add_true_lagged_effects(scen_rows, Phi = Phi)
  }

  do.call(rbind, out_list)
}

# helper: names of objects that the parallel runner needs to export
# purpose:
# keeping the export list in one place makes the parallel script easier to read.
simulation_engine_exports <- function() {
  c(
    "sample_delta_1",
    "generate_Delta_constant",
    "generate_Delta_stepwise",
    "generate_Delta_stepwise_mixture",
    "simulate_panel_data",
    "build_clpm",
    "build_riclpm",
    "build_dpm",
    "residualise_panel_linearC",
    "build_xgb_confounder_matrix",
    "tune_residualise_panel_xgb",
    "residualise_panel_xgb",
    "%||%",
    "is_whole_number_scalar",
    "safe_lavinspect",
    "safe_parameter_estimates_subset",
    "validate_common_fit_inputs",
    "validate_clpm_inputs",
    "validate_free_loading_input",
    "get_used_confounders",
    "build_preprocess_info",
    "rename_columns_for_lavaan",
    "prepare_analysis_data",
    "build_model_syntax",
    "run_lavaan_fit",
    "detect_improper_solution",
    "extract_fit_diagnostics",
    "build_target_paths",
    "empty_parameter_table",
    "extract_model_parameters",
    "fit_model_once",
    "draw_bootstrap_sample",
    "extract_bootstrap_draw",
    "summarise_bootstrap_parameters",
    "bootstrap_fit",
    "fit_clpm",
    "fit_riclpm",
    "fit_dpm",
    "fit_model_set",
    "format_fit_for_simulation",
    "add_true_lagged_effects",
    "run_one_replication"
  )
}
