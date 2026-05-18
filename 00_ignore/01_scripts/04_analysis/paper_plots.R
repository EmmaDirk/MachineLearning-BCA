# ============================================================
# libraries
# ============================================================

library(here)
library(tidyverse)
library(scales)
library(viridis)
library(knitr)


# ============================================================
# load combined data only
# ============================================================

dat_all <- readRDS(
  here("02_data", "02_server_runs", "07_run", "reduced_all_conditions_combined.rds")
)


# ============================================================
# scenario lists from combined dataframe
# ============================================================

make_scenario_list_from_combined <- function(
    dat_all,
    scenario_value,
    n_order = c(2000, 1000, 300)
) {
  this_dat <- dat_all %>%
    dplyr::mutate(
      scenario_id = as.integer(scenario_id),
      N = as.integer(N)
    ) %>%
    dplyr::filter(scenario_id == scenario_value)

  if (nrow(this_dat) == 0) {
    stop("No rows found for scenario_id = ", scenario_value, ".")
  }

  missing_n <- setdiff(n_order, sort(unique(this_dat$N)))

  if (length(missing_n) > 0) {
    stop(
      "Scenario ", scenario_value,
      " is missing these N values: ",
      paste(missing_n, collapse = ", ")
    )
  }

  purrr::map(
    n_order,
    ~ this_dat %>% dplyr::filter(N == .x)
  )
}

# Scenario 4 in the paper corresponds to scenario_id == 5 in the data.
dat_s1 <- make_scenario_list_from_combined(dat_all, scenario_value = 1)
dat_s2 <- make_scenario_list_from_combined(dat_all, scenario_value = 2)
dat_s3 <- make_scenario_list_from_combined(dat_all, scenario_value = 3)
dat_s4 <- make_scenario_list_from_combined(dat_all, scenario_value = 5)


# ============================================================
# shared helpers
# ============================================================

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

num_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

int_or_na <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

chr_or_empty <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_file <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_max_flag <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

# Include all converged models: normal convergence, warning, and severe/improper flags.
# Exclude only flag 1.
included_analysis_flags <- c(0L, 2L, 3L)
improper_analysis_flags <- c(2L, 3L)

make_extended_minor_breaks <- function(n = 16) {
  force(n)

  function(limits, major = NULL) {
    scales::breaks_extended(n = n)(limits)
  }
}

model_label <- function(model_code, free_loadings = NA) {
  model_code <- clean_chr(model_code)

  dplyr::case_when(
    model_code == "C" & !is.na(free_loadings) & free_loadings == 1 ~ "fCLPM",
    model_code == "C" ~ "CLPM",
    model_code == "R" & !is.na(free_loadings) & free_loadings == 1 ~ "fRICLPM",
    model_code == "R" ~ "RICLPM",
    model_code == "D" & !is.na(free_loadings) & free_loadings == 1 ~ "fDPM",
    model_code == "D" ~ "DPM",
    TRUE ~ model_code
  )
}

resid_label <- function(resid_code) {
  resid_code <- clean_chr(resid_code)

  dplyr::case_when(
    resid_code == "N" ~ "None",
    resid_code == "L" ~ "LM",
    resid_code == "E" ~ "EN",
    resid_code == "X" ~ "BCA",
    TRUE ~ resid_code
  )
}

direction_label <- function(effect_name) {
  dplyr::case_when(
    effect_name == "ARXY" ~ "X to Y",
    effect_name == "ARYX" ~ "Y to X",
    effect_name == "ARX"  ~ "X",
    effect_name == "ARY"  ~ "Y",
    TRUE ~ effect_name
  )
}

true_col_map <- c(
  ARX = "beta_x",
  ARY = "beta_y",
  CXY = "gamma_xy",
  CYX = "gamma_yx"
)


# ============================================================
# default model selection and method order
# ============================================================

default_keep_model_names <- c(
  "clpm_linear_confounders",
  "clpm_xgb_residualized",
  "riclpm_no_adjustment",
  "riclpm_linear_confounders",
  "riclpm_xgb_residualized"
)

default_method_order <- c(
  "CLPM adj",
  "CLPM BCA",
  "RICLPM",
  "RICLPM adj",
  "RICLPM BCA"
)


# ============================================================
# explicit method naming helper
# ============================================================

make_explicit_method_labels <- function(df) {
  df %>%
    dplyr::mutate(
      method_base = dplyr::case_when(
        model_name == "clpm_no_adjustment" ~
          "CLPM",

        model_name == "clpm_linear_confounders" ~
          "CLPM adj",

        model_name == "clpm_xgb_residualized" ~
          "CLPM BCA",

        model_name == "clpm_linear_residualized" ~
          "CLPM linear residuals",

        model_name == "clpm_enet_residualized" ~
          "CLPM EN residuals",

        model_name == "riclpm_no_adjustment" ~
          "RICLPM",

        model_name == "riclpm_linear_confounders" ~
          "RICLPM adj",

        model_name == "riclpm_xgb_residualized" ~
          "RICLPM BCA",

        model_name == "riclpm_linear_residualized" ~
          "RICLPM linear residuals",

        model_name == "riclpm_enet_residualized" ~
          "RICLPM EN residuals",

        model_name == "dpm_no_adjustment" ~
          "DPM",

        model_name == "dpm_linear_confounders" ~
          "DPM adjusted ADE",

        model_name == "dpm_xgb_residualized" ~
          "DPM BCA",

        model_name == "dpm_linear_residualized" ~
          "DPM linear residuals",

        model_name == "dpm_enet_residualized" ~
          "DPM EN residuals",

        residualizer == "N" & has_sem_adjustment ~
          paste0(model_tag, " adjusted ADE"),

        residualizer == "N" & !has_sem_adjustment ~
          model_tag,

        residualizer == "X" ~
          paste0(model_tag, " BCA"),

        residualizer == "L" ~
          paste0(model_tag, " linear residuals"),

        residualizer == "E" ~
          paste0(model_tag, " EN residuals"),

        TRUE ~ model_name
      ),

      method_detail = paste0(
        method_base,
        " [",
        dplyr::case_when(
          has_sem_adjustment ~ sem_adjustment_detail,
          has_residualizer_adjustment ~ residualizer_adjustment_detail,
          TRUE ~ "no C adjustment"
        ),
        "; ",
        se_detail,
        "]"
      ),

      method = method_base
    ) %>%
    dplyr::select(-method_base)
}


# ============================================================
# prepare data for plots and flag tables
# ============================================================

prepare_performance_data <- function(
    dat_list,
    effects_needed,
    effect_map,
    n_order = c(2000, 1000, 300),
    n_labels = c(
      `2000` = "N = 2000",
      `1000` = "N = 1000",
      `300`  = "N = 300"
    ),
    occasions = NULL,
    palette = NULL,
    keep_model_names = default_keep_model_names,
    method_order = default_method_order
) {

  stopifnot(is.list(dat_list))
  stopifnot(all(purrr::map_lgl(dat_list, is.data.frame)))

  df <- dplyr::bind_rows(dat_list) %>%
    dplyr::mutate(.source_order = dplyr::row_number())

  if (!"bootstrap_B" %in% names(df)) {
    df$bootstrap_B <- NA_integer_
  }

  missing_effects <- setdiff(effects_needed, names(effect_map))
  if (length(missing_effects) > 0) {
    stop(
      "These effects are not present in effect_map: ",
      paste(missing_effects, collapse = ", ")
    )
  }

  estimate_cols <- unique(unname(effect_map[effects_needed]))
  se_cols <- paste0("se_", estimate_cols)

  missing_true_map <- setdiff(estimate_cols, names(true_col_map))
  if (length(missing_true_map) > 0) {
    stop(
      "These estimate columns are not present in true_col_map: ",
      paste(missing_true_map, collapse = ", ")
    )
  }

  true_cols <- unique(unname(true_col_map[estimate_cols]))

  required_cols <- unique(c(
    "R", "T", "N", "analysis_flag",
    "model_name",
    "model", "residualizer",
    "sem_exclusion", "sem_c_order",
    "residualizer_exclusion", "residualizer_c_order",
    "free_loadings",
    "bootstrap_B",
    estimate_cols,
    se_cols,
    true_cols
  ))

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.null(occasions)) {
    df <- df %>%
      dplyr::filter(T %in% occasions)
  }

  oof_cols <- intersect(c("mse_x", "mse_y", "r2_x", "r2_y"), names(df))

  if (length(oof_cols) > 0) {
    df$has_oof_row <- rowSums(!is.na(df[, oof_cols, drop = FALSE])) > 0
  } else {
    df$has_oof_row <- FALSE
  }

  df <- df %>%
    dplyr::mutate(
      analysis_flag = int_or_na(analysis_flag),

      model_name = clean_chr(model_name),
      model = clean_chr(model),
      residualizer = clean_chr(residualizer),
      sem_exclusion = clean_chr(sem_exclusion),
      residualizer_exclusion = clean_chr(residualizer_exclusion),

      model_tag = model_label(model, free_loadings),
      residualizer_tag = resid_label(residualizer),

      family = dplyr::case_when(
        model == "C" ~ "CLPM",
        model == "R" ~ "RICLPM",
        model == "D" ~ "DPM",
        TRUE ~ "Other"
      ),

      sem_c_num = num_or_na(sem_c_order),
      resid_c_num = num_or_na(residualizer_c_order),
      bootstrap_B_num = num_or_na(bootstrap_B),

      has_sem_adjustment =
        residualizer == "N" &
        (
          (!is.na(sem_c_num) & sem_c_num > 0) |
            (!is.na(sem_exclusion) & sem_exclusion != "")
        ),

      has_residualizer_adjustment =
        residualizer != "N" &
        (
          (!is.na(resid_c_num) & resid_c_num > 0) |
            (!is.na(residualizer_exclusion) & residualizer_exclusion != "")
        ),

      has_confounder_adjustment =
        has_sem_adjustment | has_residualizer_adjustment,

      uses_bootstrap_se =
        !is.na(bootstrap_B_num) & bootstrap_B_num > 0,

      has_oof_row2 =
        has_oof_row |
          stringr::str_detect(
            chr_or_empty(model_name),
            stringr::regex("oof|out.?of.?fold|cross.?fit|crossfit", ignore_case = TRUE)
          ),

      sem_adjustment_detail = dplyr::case_when(
        has_sem_adjustment & is.na(sem_exclusion) ~
          paste0("SEM C order ", sem_c_num),
        has_sem_adjustment & !is.na(sem_exclusion) ~
          paste0("SEM C order ", sem_c_num, ", excluding ", sem_exclusion),
        TRUE ~ NA_character_
      ),

      residualizer_adjustment_detail = dplyr::case_when(
        has_residualizer_adjustment & is.na(residualizer_exclusion) ~
          paste0(residualizer_tag, " C order ", resid_c_num),
        has_residualizer_adjustment & !is.na(residualizer_exclusion) ~
          paste0(residualizer_tag, " C order ", resid_c_num, ", excluding ", residualizer_exclusion),
        TRUE ~ NA_character_
      ),

      se_detail = dplyr::case_when(
        uses_bootstrap_se ~ paste0("bootstrap SE, B = ", bootstrap_B_num),
        TRUE ~ "model-based SE"
      )
    ) %>%
    make_explicit_method_labels()

  if (!is.null(keep_model_names)) {
    df <- df %>%
      dplyr::filter(model_name %in% keep_model_names)
  }

  if (nrow(df) == 0) {
    stop("No rows remain after filtering by keep_model_names.")
  }

  method_id_cols <- c(
    "model_name",
    "model",
    "residualizer",
    "sem_exclusion",
    "sem_c_order",
    "residualizer_exclusion",
    "residualizer_c_order",
    "free_loadings",
    "bootstrap_B"
  )

  df <- df %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(method_id_cols), as.character)
    ) %>%
    tidyr::unite(
      col = "method_id",
      dplyr::all_of(method_id_cols),
      sep = " | ",
      remove = FALSE,
      na.rm = FALSE
    )

  method_key <- df %>%
    dplyr::group_by(
      method_id,
      method,
      method_detail,
      model_name,
      model,
      model_tag,
      residualizer,
      residualizer_tag,
      sem_exclusion,
      sem_c_order,
      residualizer_exclusion,
      residualizer_c_order,
      free_loadings,
      bootstrap_B,
      family
    ) %>%
    dplyr::summarise(
      first_row = min(.source_order, na.rm = TRUE),
      has_any_oof = any(has_oof_row2, na.rm = TRUE),
      uses_bootstrap_se = any(uses_bootstrap_se, na.rm = TRUE),
      has_sem_adjustment = any(has_sem_adjustment, na.rm = TRUE),
      has_residualizer_adjustment = any(has_residualizer_adjustment, na.rm = TRUE),
      has_confounder_adjustment = any(has_confounder_adjustment, na.rm = TRUE),
      n_rows = dplyr::n(),
      n_replications = dplyr::n_distinct(R),
      N_values = paste(sort(unique(N)), collapse = ", "),
      .groups = "drop"
    )

  method_levels <- c(
    method_order[method_order %in% unique(method_key$method)],
    setdiff(unique(method_key$method), method_order)
  )

  df <- df %>%
    dplyr::mutate(
      N_key = as.character(N),
      N_label = dplyr::coalesce(
        unname(n_labels[N_key]),
        paste0("N = ", N_key)
      )
    )

  n_levels <- unname(n_labels[as.character(n_order)])
  missing_n_levels <- is.na(n_levels)
  n_levels[missing_n_levels] <- paste0("N = ", n_order[missing_n_levels])

  df <- df %>%
    dplyr::mutate(
      method = factor(method, levels = method_levels),
      N_label = factor(N_label, levels = n_levels)
    )

  if (is.null(palette)) {
    if (requireNamespace("viridis", quietly = TRUE)) {
      palette <- setNames(
        viridis::viridis(
          n = length(method_levels),
          option = "C",
          end = 0.9
        ),
        method_levels
      )
    } else {
      palette <- setNames(
        grDevices::hcl.colors(
          n = length(method_levels),
          palette = "viridis"
        ),
        method_levels
      )
    }
  } else {
    missing_palette_names <- setdiff(method_levels, names(palette))

    if (length(missing_palette_names) > 0) {
      if (requireNamespace("viridis", quietly = TRUE)) {
        extra_cols <- viridis::viridis(
          n = length(missing_palette_names),
          option = "C",
          end = 0.9
        )
      } else {
        extra_cols <- grDevices::hcl.colors(
          n = length(missing_palette_names),
          palette = "viridis"
        )
      }

      palette <- c(palette, setNames(extra_cols, missing_palette_names))
    }

    palette <- palette[method_levels]
  }

  model_overview <- method_key %>%
    dplyr::select(
      method,
      method_detail,
      model_name,
      model,
      residualizer,
      residualizer_tag,
      sem_c_order,
      residualizer_c_order,
      sem_exclusion,
      residualizer_exclusion,
      free_loadings,
      bootstrap_B,
      first_row,
      has_any_oof,
      uses_bootstrap_se,
      has_sem_adjustment,
      has_residualizer_adjustment,
      has_confounder_adjustment,
      n_replications,
      N_values
    ) %>%
    dplyr::arrange(factor(method, levels = method_levels), model_name, first_row)

  list(
    df = df,
    method_key = method_key,
    method_levels = method_levels,
    n_levels = n_levels,
    palette = palette,
    model_overview = model_overview
  )
}


# ============================================================
# plot function
# ============================================================

plot_performance_effects_3x4 <- function(
    dat_list,
    effects = c("ARXY"),
    effect_map = c(
      ARX  = "ARX",
      ARY  = "ARY",
      ARXY = "CXY",
      ARYX = "CYX"
    ),
    type1_effect_for = c(
      ARXY = "ARYX"
    ),
    type1_metric_labels = c(
      ARXY = "Type-I error Y to X"
    ),
    # Changed for plot facets:
    # columns now appear as N = 300, N = 1000, N = 2000.
    n_order = c(300, 1000, 2000),
    n_labels = c(
      `300`  = "N = 300",
      `1000` = "N = 1000",
      `2000` = "N = 2000"
    ),
    occasions = 2:5,
    alpha = 0.05,
    drop_flagged = TRUE,
    width = 12,
    height = 16.5,
    y_axis_breaks = 5,
    y_axis_minor_breaks = 16,
    include_mcse_bars = TRUE,
    palette = NULL,
    keep_model_names = default_keep_model_names,
    method_order = default_method_order,
    use_relative_bias_when_true_nonzero = TRUE
) {

  type1_effect_for <- type1_effect_for[!is.na(type1_effect_for)]
  type1_effect_for <- type1_effect_for[names(type1_effect_for) %in% effects]

  effects_needed <- unique(c(effects, unname(type1_effect_for)))

  prep <- prepare_performance_data(
    dat_list = dat_list,
    effects_needed = effects_needed,
    effect_map = effect_map,
    n_order = n_order,
    n_labels = n_labels,
    occasions = occasions,
    palette = palette,
    keep_model_names = keep_model_names,
    method_order = method_order
  )

  df <- prep$df
  method_levels <- prep$method_levels
  n_levels <- prep$n_levels
  palette <- prep$palette

  keep_valid_estimates <- function(x) {
    if (!drop_flagged) {
      return(x)
    }

    x %>%
      dplyr::filter(analysis_flag %in% included_analysis_flags)
  }

  make_sim_long_from_df <- function(this_effect) {
    estimate_col <- unname(effect_map[this_effect])
    se_col <- paste0("se_", estimate_col)
    true_col <- unname(true_col_map[estimate_col])

    df %>%
      dplyr::transmute(
        effect = this_effect,
        R,
        T,
        N,
        N_label,
        analysis_flag,
        model_name,
        method_id,
        method,
        method_detail,
        family,
        estimate = .data[[estimate_col]],
        se_estimate = .data[[se_col]],
        true = .data[[true_col]],
        error = estimate - true,
        abs_error = abs(error),
        sq_error = error^2
      )
  }

  build_one_effect <- function(effect_name) {

    main_direction <- direction_label(effect_name)

    sim_long <- make_sim_long_from_df(effect_name)
    sim_long <- keep_valid_estimates(sim_long)

    true_values <- sim_long %>%
      dplyr::filter(!is.na(true)) %>%
      dplyr::pull(true)

    uses_raw_bias <-
      !use_relative_bias_when_true_nonzero ||
      length(true_values) == 0 ||
      any(abs(true_values) < 1e-12)

    bias_label <- if (uses_raw_bias) {
      "Bias"
    } else {
      "Relative Bias"
    }

    extra_type1_effect <- NA_character_
    extra_type1_label <- NULL

    if (!is.null(type1_effect_for) && effect_name %in% names(type1_effect_for)) {
      extra_type1_effect <- unname(type1_effect_for[effect_name])

      extra_type1_label <- if (
        !is.null(type1_metric_labels) &&
          effect_name %in% names(type1_metric_labels)
      ) {
        unname(type1_metric_labels[effect_name])
      } else {
        paste0("Type-I error ", direction_label(extra_type1_effect))
      }
    }

    metric_labels <- c(
      bias = paste(bias_label, main_direction),
      rmse = paste("RMSE", main_direction),
      power = paste("Power", main_direction)
    )

    if (!is.null(extra_type1_label)) {
      metric_labels <- c(
        metric_labels,
        type1_extra = unname(extra_type1_label)
      )
    }

    if (uses_raw_bias) {
      bias_df <- sim_long %>%
        dplyr::filter(!is.na(estimate), !is.na(true)) %>%
        dplyr::group_by(N, N_label, T, method_id, method, method_detail, family) %>%
        dplyr::summarise(
          nsim = dplyr::n(),
          true = dplyr::first(true),
          mean_est = mean(estimate, na.rm = TRUE),
          value = mean_est - true,
          mcse = stats::sd(estimate, na.rm = TRUE) / sqrt(nsim),
          .groups = "drop"
        ) %>%
        dplyr::mutate(metric = "bias")
    } else {
      bias_df <- sim_long %>%
        dplyr::filter(!is.na(estimate), !is.na(true), abs(true) > 1e-12) %>%
        dplyr::group_by(N, N_label, T, method_id, method, method_detail, family) %>%
        dplyr::summarise(
          nsim = dplyr::n(),
          true = dplyr::first(true),
          mean_est = mean(estimate, na.rm = TRUE),
          value = (mean_est - true) / true,
          mcse = stats::sd(estimate, na.rm = TRUE) / sqrt(nsim) / abs(true),
          .groups = "drop"
        ) %>%
        dplyr::mutate(metric = "bias")
    }

    rmse_df <- sim_long %>%
      dplyr::filter(!is.na(estimate), !is.na(true)) %>%
      dplyr::mutate(sq_err = (estimate - true)^2) %>%
      dplyr::group_by(N, N_label, T, method_id, method, method_detail, family) %>%
      dplyr::summarise(
        nsim = dplyr::n(),
        mean_sq_err = mean(sq_err, na.rm = TRUE),
        var_sq_err = stats::var(sq_err, na.rm = TRUE),
        value = sqrt(mean_sq_err),
        mcse = dplyr::if_else(
          nsim <= 1 | is.na(value) | value == 0,
          NA_real_,
          sqrt(var_sq_err / nsim) / (2 * value)
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(metric = "rmse")

    crit_z <- stats::qnorm(1 - alpha / 2)

    make_detection_df <- function(this_sim_long, metric_name) {
      this_sim_long %>%
        dplyr::filter(
          !is.na(estimate),
          !is.na(se_estimate),
          se_estimate > 0,
          !is.na(true)
        ) %>%
        dplyr::mutate(
          z_value = estimate / se_estimate,
          detected = abs(z_value) > crit_z
        ) %>%
        dplyr::group_by(N, N_label, T, method_id, method, method_detail, family) %>%
        dplyr::summarise(
          nsim = dplyr::n(),
          true = dplyr::first(true),
          detect_prob = mean(detected, na.rm = TRUE),
          value = detect_prob,
          mcse = sqrt(detect_prob * (1 - detect_prob) / nsim),
          detection_type = dplyr::if_else(abs(true) < 1e-12, "Type-I error", "Power"),
          .groups = "drop"
        ) %>%
        dplyr::mutate(metric = metric_name)
    }

    power_df <- make_detection_df(
      this_sim_long = sim_long,
      metric_name = "power"
    )

    extra_type1_df <- NULL

    if (!is.na(extra_type1_effect)) {
      sim_long_type1 <- make_sim_long_from_df(extra_type1_effect)
      sim_long_type1 <- keep_valid_estimates(sim_long_type1)

      extra_type1_df <- make_detection_df(
        this_sim_long = sim_long_type1,
        metric_name = "type1_extra"
      )
    }

    plot_df <- dplyr::bind_rows(
      bias_df,
      rmse_df,
      power_df,
      extra_type1_df
    ) %>%
      dplyr::mutate(
        N_label = factor(N_label, levels = n_levels),
        method = factor(method, levels = method_levels),

        ymin = dplyr::case_when(
          metric %in% c("rmse", "power", "type1_extra") ~ pmax(0, value - mcse),
          TRUE ~ value - mcse
        ),

        ymax = dplyr::case_when(
          metric %in% c("power", "type1_extra") ~ pmin(1, value + mcse),
          TRUE ~ value + mcse
        ),

        metric = factor(
          metric,
          levels = names(metric_labels),
          labels = unname(metric_labels)
        )
      )

    ref_df <- tibble::tibble(
      metric = "bias",
      ref_value = 0
    )

    if (!is.null(extra_type1_label)) {
      ref_df <- dplyr::bind_rows(
        ref_df,
        tibble::tibble(
          metric = "type1_extra",
          ref_value = alpha
        )
      )
    }

    ref_df <- ref_df %>%
      dplyr::mutate(
        metric = factor(
          metric,
          levels = names(metric_labels),
          labels = unname(metric_labels)
        )
      )

    panel_bg_df <- plot_df %>%
      dplyr::distinct(metric, N_label) %>%
      dplyr::mutate(
        metric_index = as.integer(metric),
        panel_fill = dplyr::if_else(metric_index %% 2 == 1, "white", "grey96")
      )

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = T,
        y = value,
        color = method,
        group = method_id
      )
    ) +
      ggplot2::geom_rect(
        data = panel_bg_df,
        ggplot2::aes(
          xmin = -Inf,
          xmax = Inf,
          ymin = -Inf,
          ymax = Inf,
          fill = panel_fill
        ),
        inherit.aes = FALSE,
        color = NA
      ) +
      ggplot2::scale_fill_identity(guide = "none") +
      ggplot2::geom_hline(
        data = ref_df,
        ggplot2::aes(yintercept = ref_value),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.5
      ) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::geom_point(size = 2.7)

    if (include_mcse_bars) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ymin, ymax = ymax),
          width = 0.15,
          linewidth = 0.55
        )
    }

    p <- p +
      ggplot2::facet_grid(metric ~ N_label, scales = "free_y") +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_extended(n = y_axis_breaks),
        minor_breaks = make_extended_minor_breaks(n = y_axis_minor_breaks)
      ) +
      ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$T))) +
      ggplot2::scale_color_manual(
        values = palette,
        limits = method_levels,
        breaks = method_levels,
        drop = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          nrow = 1,
          byrow = TRUE,
          override.aes = list(linewidth = 1.2, size = 3.2)
        )
      ) +
      ggplot2::labs(
        title = NULL,
        x = "Occasion",
        y = NULL,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),

        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.border = ggplot2::element_rect(fill = NA, color = "grey80", linewidth = 0.4),
        panel.spacing = grid::unit(1.15, "lines"),

        panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor.y = ggplot2::element_line(color = "grey92", linewidth = 0.2),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),

        strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
        strip.text = ggplot2::element_text(size = 15, face = "plain"),

        axis.title = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(size = 13),

        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 13),
        legend.key.width = grid::unit(1.8, "lines"),
        legend.key.height = grid::unit(1.05, "lines"),
        legend.spacing.y = grid::unit(0.25, "lines"),
        legend.margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0),
        legend.box = "horizontal"
      )

    list(
      plot = p,
      plot_df = plot_df,
      sim_long = sim_long
    )
  }

  effect_results <- purrr::map(effects, build_one_effect)
  names(effect_results) <- effects

  plots <- purrr::map(effect_results, "plot")
  plot_dfs <- purrr::map(effect_results, "plot_df")

  metric_df <- dplyr::bind_rows(plot_dfs, .id = "requested_effect")

  analysis_long <- purrr::map_dfr(effects_needed, make_sim_long_from_df)

  analysis_long_used <- if (drop_flagged) {
    analysis_long %>%
      dplyr::filter(analysis_flag %in% included_analysis_flags)
  } else {
    analysis_long
  }

  list(
    plots = plots,
    plot_dfs = plot_dfs,
    metric_df = metric_df,
    raw_df = df,
    analysis_long = analysis_long,
    analysis_long_used = analysis_long_used,
    effect_results = effect_results,
    model_overview = prep$model_overview,
    method_key = prep$method_key,
    method_levels = method_levels,
    n_levels = n_levels,
    palette = palette
  )
}


# ============================================================
# improper fit proportions and combined latex table
# ============================================================

make_improper_fit_table <- function(
    dat_list,
    n_order = c(2000, 1000, 300),
    keep_model_names = default_keep_model_names,
    method_order = default_method_order
) {
  prep <- prepare_performance_data(
    dat_list = dat_list,
    effects_needed = c("ARXY", "ARYX"),
    effect_map = c(
      ARX  = "ARX",
      ARY  = "ARY",
      ARXY = "CXY",
      ARYX = "CYX"
    ),
    n_order = n_order,
    occasions = NULL,
    keep_model_names = keep_model_names,
    method_order = method_order
  )

  prep$df %>%
    dplyr::mutate(analysis_flag = int_or_na(analysis_flag)) %>%
    dplyr::group_by(N, method, method_id, model_name, R) %>%
    dplyr::summarise(
      original_fit_code = safe_max_flag(analysis_flag),
      .groups = "drop"
    ) %>%
    dplyr::group_by(N, method) %>%
    dplyr::summarise(
      n_replications = dplyr::n_distinct(R),
      improper_n = sum(original_fit_code %in% improper_analysis_flags, na.rm = TRUE),
      improper_frac = improper_n / n_replications,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      N = as.integer(N),
      method = factor(method, levels = method_order)
    ) %>%
    dplyr::arrange(match(N, n_order), method)
}

format_latex_number <- function(x) {
  sprintf("%.2f", x)
}

make_latex_improper_table_all_scenarios <- function(
    improper_list,
    scenario_numbers = names(improper_list),
    n_order = c(2000, 1000, 300),
    method_cols = default_method_order,
    caption = "Fraction of Improper Solutions by Sample Size, Method, and Scenario",
    label = "tab:improper-all-scenarios"
) {
  if (is.null(scenario_numbers) || any(scenario_numbers == "")) {
    scenario_numbers <- seq_along(improper_list)
  }

  combined_df <- purrr::map2_dfr(
    improper_list,
    scenario_numbers,
    ~ .x %>%
      dplyr::mutate(Scenario = as.character(.y))
  )

  wide_df <- combined_df %>%
    dplyr::mutate(method = as.character(method)) %>%
    dplyr::select(Scenario, N, method, improper_frac) %>%
    tidyr::pivot_wider(
      names_from = method,
      values_from = improper_frac
    )

  missing_methods <- setdiff(method_cols, names(wide_df))
  if (length(missing_methods) > 0) {
    wide_df[missing_methods] <- NA_real_
  }

  table_df <- wide_df %>%
    dplyr::select(Scenario, N, dplyr::all_of(method_cols)) %>%
    dplyr::mutate(
      Scenario = factor(Scenario, levels = as.character(scenario_numbers)),
      N = as.integer(N)
    ) %>%
    dplyr::arrange(Scenario, match(N, n_order)) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(method_cols),
        ~ dplyr::if_else(is.na(.x), "--", format_latex_number(.x))
      )
    )

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{ll", paste(rep("c", length(method_cols)), collapse = ""), "}"),
    "\\toprule",
    paste0("Scenario & $N$ & ", paste(method_cols, collapse = " & "), "\\\\"),
    "\\midrule"
  )

  scenario_levels <- as.character(scenario_numbers)

  for (i in seq_along(scenario_levels)) {
    this_scenario <- scenario_levels[i]

    this_df <- table_df %>%
      dplyr::filter(as.character(Scenario) == this_scenario) %>%
      dplyr::arrange(match(N, n_order))

    for (j in seq_len(nrow(this_df))) {
      scenario_cell <- if (j == 1) this_scenario else " "
      n_cell <- as.character(this_df$N[j])
      value_cells <- unname(unlist(this_df[j, method_cols], use.names = FALSE))

      lines <- c(
        lines,
        paste0(
          paste(c(scenario_cell, n_cell, value_cells), collapse = " & "),
          "\\\\"
        )
      )
    }

    if (i < length(scenario_levels)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  paste(lines, collapse = "\n")
}


# ============================================================
# create plots
# ============================================================

perf_s1 <- plot_performance_effects_3x4(dat_list = dat_s1)
perf_s2 <- plot_performance_effects_3x4(dat_list = dat_s2)
perf_s3 <- plot_performance_effects_3x4(dat_list = dat_s3)
perf_s4 <- plot_performance_effects_3x4(dat_list = dat_s4)


# ============================================================
# export final figures only
# ============================================================

output_dir <- here("04_plots", "final_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggplot2::ggsave(
  filename = file.path(output_dir, "scenario_1.png"),
  plot = perf_s1$plots$ARXY,
  width = 12,
  height = 16.5,
  dpi = 300,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "scenario_2.png"),
  plot = perf_s2$plots$ARXY,
  width = 12,
  height = 16.5,
  dpi = 300,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "scenario_3.png"),
  plot = perf_s3$plots$ARXY,
  width = 12,
  height = 16.5,
  dpi = 300,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "scenario_4.png"),
  plot = perf_s4$plots$ARXY,
  width = 12,
  height = 16.5,
  dpi = 300,
  units = "in",
  bg = "white"
)


# ============================================================
# export one combined latex table only
# ============================================================

improper_s1 <- make_improper_fit_table(dat_s1)
improper_s2 <- make_improper_fit_table(dat_s2)
improper_s3 <- make_improper_fit_table(dat_s3)
improper_s4 <- make_improper_fit_table(dat_s4)

latex_table_all_scenarios <- make_latex_improper_table_all_scenarios(
  improper_list = list(
    `1` = improper_s1,
    `2` = improper_s2,
    `3` = improper_s3,
    `4` = improper_s4
  )
)

writeLines(
  latex_table_all_scenarios,
  con = file.path(output_dir, "improper_fit_proportions.tex")
)


# ============================================================
# run from console, not from inside this script
# ============================================================

# source(
#   here::here("01_scripts", "04_analysis", "paper_plots.R"),
#   echo = TRUE
# )