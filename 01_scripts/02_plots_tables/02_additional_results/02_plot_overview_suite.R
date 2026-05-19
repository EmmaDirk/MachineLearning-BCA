# =================================================================================================
#
# This script contains helper functions for inspecting combined simulation-results data frames.
#
# The model-overview helpers create compact tables of the model specifications in the results file.
# The overview-suite function creates diagnostic summary tables and plots for SEM performance,
# residualizer diagnostics, true delta trajectories, convergence flags, and improper solutions.
#
# Main functions:
#   overview_models()
#   plot_overview_suite()
# =================================================================================================

# ---- small shared helpers ------------------------------------------------------------------------

sim_clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

sim_num_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

sim_int_or_na <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

sim_chr_or_empty <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

sim_fmt_excl <- function(x) {
  x <- sim_clean_chr(x)
  ifelse(is.na(x), "no", x)
}

sim_fmt_order <- function(x) {
  ifelse(is.na(x), NA_character_, as.character(x))
}

sim_model_label <- function(model_code, free_loadings = NA) {
  model_code <- sim_clean_chr(model_code)

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

sim_resid_label <- function(resid_code, x_label = "BCA") {
  resid_code <- sim_clean_chr(resid_code)

  dplyr::case_when(
    resid_code == "N" ~ "None",
    resid_code == "L" ~ "LM",
    resid_code == "E" ~ "EN",
    resid_code == "X" ~ x_label,
    TRUE ~ resid_code
  )
}

sim_direction_label <- function(effect_name) {
  dplyr::case_when(
    effect_name == "ARXY" ~ "X to Y",
    effect_name == "ARYX" ~ "Y to X",
    effect_name == "CXY"  ~ "X to Y",
    effect_name == "CYX"  ~ "Y to X",
    effect_name == "ARX"  ~ "X autoregression",
    effect_name == "ARY"  ~ "Y autoregression",
    TRUE ~ effect_name
  )
}

sim_make_id <- function(data, cols) {
  stopifnot(is.data.frame(data))

  pieces <- lapply(cols, function(col) {
    value <- if (col %in% names(data)) {
      as.character(data[[col]])
    } else {
      rep(NA_character_, nrow(data))
    }

    value[is.na(value)] <- "NA"
    paste0(col, "=", value)
  })

  do.call(paste, c(pieces, sep = " | "))
}

sim_filter_results <- function(
    results_df,
    scenario_id = NULL,
    N = NULL,
    occasions = NULL,
    keep_model_names = NULL,
    keep_families = NULL
) {
  stopifnot(is.data.frame(results_df))

  df <- results_df

  if (!is.null(scenario_id)) {
    if (!"scenario_id" %in% names(df)) {
      stop("scenario_id was supplied, but results_df has no scenario_id column.")
    }

    df <- df %>%
      dplyr::filter(.data$scenario_id %in% .env$scenario_id)
  }

  if (!is.null(N)) {
    if (!"N" %in% names(df)) {
      stop("N was supplied, but results_df has no N column.")
    }

    df <- df %>%
      dplyr::filter(.data$N %in% .env$N)
  }

  if (!is.null(occasions)) {
    if (!"T" %in% names(df)) {
      stop("occasions was supplied, but results_df has no T column.")
    }

    df <- df %>%
      dplyr::filter(.data$T %in% .env$occasions)
  }

  if (!is.null(keep_model_names)) {
    if (!"model_name" %in% names(df)) {
      stop("keep_model_names was supplied, but results_df has no model_name column.")
    }

    df <- df %>%
      dplyr::filter(.data$model_name %in% .env$keep_model_names)
  }

  if (!is.null(keep_families)) {
    if (!"model" %in% names(df)) {
      stop("keep_families was supplied, but results_df has no model column.")
    }

    wanted_codes <- dplyr::case_when(
      keep_families %in% c("CLPM", "C") ~ "C",
      keep_families %in% c("RICLPM", "RI-CLPM", "R") ~ "R",
      keep_families %in% c("DPM", "D") ~ "D",
      TRUE ~ keep_families
    )

    df <- df %>%
      dplyr::filter(.data$model %in% .env$wanted_codes)
  }

  df
}


# ---- method parsing ------------------------------------------------------------------------------

sim_make_method_columns <- function(df) {
  stopifnot(is.data.frame(df))

  n <- nrow(df)

  optional_cols <- list(
    scenario_id = NA_integer_,
    scenario_label = NA_character_,
    N = NA_integer_,
    model_name = NA_character_,
    model = NA_character_,
    residualizer = NA_character_,
    sem_exclusion = NA_character_,
    sem_c_order = NA_real_,
    residualizer_exclusion = NA_character_,
    residualizer_c_order = NA_real_,
    free_loadings = NA_real_,
    bootstrap_B = NA_real_,
    analysis_flag = 0L,
    flag0 = NA_real_,
    flag1 = NA_real_,
    flag2 = NA_real_,
    flag3 = NA_real_,
    improper_reason = NA_character_
  )

  for (nm in names(optional_cols)) {
    if (!nm %in% names(df)) {
      df[[nm]] <- rep(optional_cols[[nm]], n)
    }
  }

  if (!"R" %in% names(df)) {
    df$R <- seq_len(n)
  }

  if (!"T" %in% names(df)) {
    stop("results_df must contain a T column.")
  }

  df <- df %>%
    dplyr::mutate(
      .source_order = dplyr::row_number(),

      scenario_id = sim_int_or_na(.data$scenario_id),
      N = sim_int_or_na(.data$N),
      analysis_flag = sim_int_or_na(.data$analysis_flag),

      model_name = sim_clean_chr(.data$model_name),
      model = sim_clean_chr(.data$model),
      residualizer = sim_clean_chr(.data$residualizer),
      sem_exclusion = sim_clean_chr(.data$sem_exclusion),
      residualizer_exclusion = sim_clean_chr(.data$residualizer_exclusion),
      improper_reason = sim_clean_chr(.data$improper_reason),

      scenario_label = dplyr::case_when(
        !is.na(sim_clean_chr(.data$scenario_label)) ~ sim_clean_chr(.data$scenario_label),
        !is.na(.data$scenario_id) ~ paste0("scenario_id = ", .data$scenario_id),
        TRUE ~ NA_character_
      ),

      N_label = dplyr::case_when(
        !is.na(.data$N) ~ paste0("N = ", .data$N),
        TRUE ~ NA_character_
      ),

      model_tag = sim_model_label(.data$model, .data$free_loadings),
      residualizer_tag = sim_resid_label(.data$residualizer),

      family = dplyr::case_when(
        .data$model == "C" ~ "CLPM",
        .data$model == "R" ~ "RICLPM",
        .data$model == "D" ~ "DPM",
        TRUE ~ "Other"
      ),

      sem_c_num = sim_num_or_na(.data$sem_c_order),
      resid_c_num = sim_num_or_na(.data$residualizer_c_order),
      bootstrap_B_num = sim_num_or_na(.data$bootstrap_B),

      has_sem_adjustment =
        .data$residualizer == "N" &
        (
          (!is.na(.data$sem_c_num) & .data$sem_c_num > 0) |
            (!is.na(.data$sem_exclusion) & .data$sem_exclusion != "")
        ),

      has_residualizer_adjustment =
        .data$residualizer != "N" &
        (
          (!is.na(.data$resid_c_num) & .data$resid_c_num > 0) |
            (!is.na(.data$residualizer_exclusion) & .data$residualizer_exclusion != "")
        ),

      has_confounder_adjustment =
        .data$has_sem_adjustment | .data$has_residualizer_adjustment,

      uses_bootstrap_se =
        !is.na(.data$bootstrap_B_num) & .data$bootstrap_B_num > 0,

      sem_adjustment_detail = dplyr::case_when(
        .data$has_sem_adjustment & is.na(.data$sem_exclusion) ~
          paste0("SEM C order ", .data$sem_c_num),
        .data$has_sem_adjustment & !is.na(.data$sem_exclusion) ~
          paste0("SEM C order ", .data$sem_c_num, ", excluding ", .data$sem_exclusion),
        TRUE ~ NA_character_
      ),

      residualizer_adjustment_detail = dplyr::case_when(
        .data$has_residualizer_adjustment & is.na(.data$residualizer_exclusion) ~
          paste0(.data$residualizer_tag, " C order ", .data$resid_c_num),
        .data$has_residualizer_adjustment & !is.na(.data$residualizer_exclusion) ~
          paste0(.data$residualizer_tag, " C order ", .data$resid_c_num, ", excluding ", .data$residualizer_exclusion),
        TRUE ~ NA_character_
      ),

      se_detail = dplyr::case_when(
        .data$uses_bootstrap_se ~ paste0("bootstrap SE, B = ", .data$bootstrap_B_num),
        TRUE ~ "model-based SE"
      ),

      method_base = dplyr::case_when(
        .data$model_name == "clpm_no_adjustment" ~ "CLPM",
        .data$model_name == "clpm_linear_confounders" ~ "CLPM adj",
        .data$model_name == "clpm_xgb_residualized" ~ "CLPM BCA",
        .data$model_name == "clpm_linear_residualized" ~ "CLPM linear residuals",
        .data$model_name == "clpm_enet_residualized" ~ "CLPM EN residuals",

        .data$model_name == "riclpm_no_adjustment" ~ "RICLPM",
        .data$model_name == "riclpm_linear_confounders" ~ "RICLPM adj",
        .data$model_name == "riclpm_xgb_residualized" ~ "RICLPM BCA",
        .data$model_name == "riclpm_linear_residualized" ~ "RICLPM linear residuals",
        .data$model_name == "riclpm_enet_residualized" ~ "RICLPM EN residuals",

        .data$model_name == "dpm_no_adjustment" ~ "DPM",
        .data$model_name == "dpm_linear_confounders" ~ "DPM adjusted ADE",
        .data$model_name == "dpm_xgb_residualized" ~ "DPM BCA",
        .data$model_name == "dpm_linear_residualized" ~ "DPM linear residuals",
        .data$model_name == "dpm_enet_residualized" ~ "DPM EN residuals",

        .data$residualizer == "N" & .data$has_sem_adjustment ~
          paste0(.data$model_tag, " adjusted ADE"),

        .data$residualizer == "N" & !.data$has_sem_adjustment ~
          .data$model_tag,

        .data$residualizer == "X" ~
          paste0(.data$model_tag, " BCA"),

        .data$residualizer == "L" ~
          paste0(.data$model_tag, " linear residuals"),

        .data$residualizer == "E" ~
          paste0(.data$model_tag, " EN residuals"),

        !is.na(.data$model_name) ~ .data$model_name,
        TRUE ~ paste0(.data$model_tag, " / ", .data$residualizer_tag)
      ),

      method_detail = paste0(
        .data$method_base,
        " [",
        dplyr::case_when(
          .data$has_sem_adjustment ~ .data$sem_adjustment_detail,
          .data$has_residualizer_adjustment ~ .data$residualizer_adjustment_detail,
          TRUE ~ "no C adjustment"
        ),
        "; ",
        .data$se_detail,
        "]"
      ),

      method = .data$method_base
    ) %>%
    dplyr::select(-method_base)

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

  df$method_id <- sim_make_id(df, method_id_cols)

  df
}


# ---- overview_models -----------------------------------------------------------------------------

overview_models <- function(
    results_df,
    scenario_id = NULL,
    N = NULL,
    occasions = NULL,
    keep_model_names = NULL,
    keep_families = NULL,
    include_counts = TRUE,
    print_result = TRUE
) {
  stopifnot(is.data.frame(results_df))

  df <- results_df %>%
    sim_filter_results(
      scenario_id = scenario_id,
      N = N,
      occasions = occasions,
      keep_model_names = keep_model_names,
      keep_families = keep_families
    ) %>%
    sim_make_method_columns()

  if (nrow(df) == 0) {
    stop("No rows remain after filtering.")
  }

  model_overview <- df %>%
    dplyr::group_by(
      .data$method_id,
      .data$method,
      .data$method_detail,
      .data$family,
      .data$model_name,
      .data$model,
      .data$model_tag,
      .data$residualizer,
      .data$residualizer_tag,
      .data$sem_exclusion,
      .data$sem_c_order,
      .data$residualizer_exclusion,
      .data$residualizer_c_order,
      .data$free_loadings,
      .data$bootstrap_B
    ) %>%
    dplyr::summarise(
      first_row = min(.data$.source_order, na.rm = TRUE),
      n_rows = dplyr::n(),
      n_replications = dplyr::n_distinct(.data$R),
      T_values = paste(sort(unique(.data$T)), collapse = ", "),
      N_values = paste(sort(unique(.data$N)), collapse = ", "),
      scenario_values = paste(sort(unique(.data$scenario_id)), collapse = ", "),
      uses_bootstrap_se = any(.data$uses_bootstrap_se, na.rm = TRUE),
      has_sem_adjustment = any(.data$has_sem_adjustment, na.rm = TRUE),
      has_residualizer_adjustment = any(.data$has_residualizer_adjustment, na.rm = TRUE),
      has_confounder_adjustment = any(.data$has_confounder_adjustment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$family, .data$method, .data$model_name, .data$first_row)

  display_cols <- c(
    "method",
    "method_detail",
    "family",
    "model_name",
    "model",
    "residualizer",
    "residualizer_tag",
    "sem_c_order",
    "residualizer_c_order",
    "sem_exclusion",
    "residualizer_exclusion",
    "free_loadings",
    "bootstrap_B",
    "uses_bootstrap_se",
    "has_sem_adjustment",
    "has_residualizer_adjustment",
    "has_confounder_adjustment"
  )

  if (include_counts) {
    display_cols <- c(
      display_cols,
      "n_rows",
      "n_replications",
      "T_values",
      "N_values",
      "scenario_values"
    )
  }

  model_overview <- model_overview %>%
    dplyr::select(dplyr::any_of(display_cols)) %>%
    tibble::as_tibble()

  if (print_result) {
    print(model_overview, n = Inf, width = Inf)
  }

  invisible(model_overview)
}


# ---- defaults ------------------------------------------------------------------------------------

default_true_col_map <- c(
  ARX = "beta_x",
  ARY = "beta_y",
  CXY = "gamma_xy",
  CYX = "gamma_yx"
)

default_effect_map <- c(
  ARX = "ARX",
  ARY = "ARY",
  CXY = "CXY",
  CYX = "CYX",
  ARXY = "CXY",
  ARYX = "CYX"
)

default_included_analysis_flags <- c(0L, 2L, 3L)


# ---- generic helpers -----------------------------------------------------------------------------

summarise_mean_mcse <- function(data, value_col, mean_name, mcse_name, group_cols) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      !!mean_name := mean(.data[[value_col]], na.rm = TRUE),
      !!mcse_name := stats::sd(.data[[value_col]], na.rm = TRUE) / sqrt(.data$nsim),
      .groups = "drop"
    )
}

infer_ordered_levels <- function(x, preferred = NULL) {
  x <- as.character(x)
  x <- x[!is.na(x)]

  c(
    preferred[preferred %in% unique(x)],
    setdiff(unique(x), preferred)
  )
}

make_palette <- function(levels, palette = NULL) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels)]

  if (length(levels) == 0) {
    return(stats::setNames(character(0), character(0)))
  }

  if (is.null(palette)) {
    palette <- stats::setNames(
      scales::hue_pal()(length(levels)),
      levels
    )
  } else {
    missing_names <- setdiff(levels, names(palette))

    if (length(missing_names) > 0) {
      extra_cols <- scales::hue_pal()(length(missing_names))
      palette <- c(palette, stats::setNames(extra_cols, missing_names))
    }

    palette <- palette[levels]
  }

  palette
}

prepare_suite_data <- function(
    results_df,
    scenario_id = NULL,
    N = NULL,
    occasions = NULL,
    keep_model_names = NULL,
    keep_families = NULL,
    method_order = NULL,
    n_order = NULL,
    n_labels = NULL,
    palette = NULL,
    print_messages = TRUE
) {
  df <- results_df %>%
    sim_filter_results(
      scenario_id = scenario_id,
      N = N,
      occasions = occasions,
      keep_model_names = keep_model_names,
      keep_families = keep_families
    ) %>%
    sim_make_method_columns()

  if (nrow(df) == 0) {
    stop("No rows remain after filtering.")
  }

  if (is.null(n_order)) {
    n_order <- sort(unique(df$N[!is.na(df$N)]))
  }

  if (is.null(n_labels)) {
    n_labels <- stats::setNames(paste0("N = ", n_order), as.character(n_order))
  }

  n_levels <- unname(n_labels[as.character(n_order)])
  missing_n_levels <- is.na(n_levels)

  if (any(missing_n_levels)) {
    n_levels[missing_n_levels] <- paste0("N = ", n_order[missing_n_levels])
  }

  df <- df %>%
    dplyr::mutate(
      N_key = as.character(.data$N),
      N_label = dplyr::coalesce(unname(n_labels[.data$N_key]), .data$N_label),
      N_label = factor(.data$N_label, levels = n_levels)
    )

  method_key <- df %>%
    dplyr::group_by(
      .data$method_id,
      .data$method,
      .data$method_detail,
      .data$model_name,
      .data$model,
      .data$model_tag,
      .data$residualizer,
      .data$residualizer_tag,
      .data$sem_exclusion,
      .data$sem_c_order,
      .data$residualizer_exclusion,
      .data$residualizer_c_order,
      .data$free_loadings,
      .data$bootstrap_B,
      .data$family
    ) %>%
    dplyr::summarise(
      first_row = min(.data$.source_order, na.rm = TRUE),
      uses_bootstrap_se = any(.data$uses_bootstrap_se, na.rm = TRUE),
      has_sem_adjustment = any(.data$has_sem_adjustment, na.rm = TRUE),
      has_residualizer_adjustment = any(.data$has_residualizer_adjustment, na.rm = TRUE),
      has_confounder_adjustment = any(.data$has_confounder_adjustment, na.rm = TRUE),
      n_rows = dplyr::n(),
      n_replications = dplyr::n_distinct(.data$R),
      N_values = paste(sort(unique(.data$N)), collapse = ", "),
      scenario_values = paste(sort(unique(.data$scenario_id)), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$family, .data$method, .data$first_row)

  method_levels <- infer_ordered_levels(method_key$method, method_order)
  pal_method <- make_palette(method_levels, palette)

  df <- df %>%
    dplyr::mutate(
      method = factor(as.character(.data$method), levels = method_levels)
    )

  method_key <- method_key %>%
    dplyr::mutate(
      method = factor(as.character(.data$method), levels = method_levels)
    )

  model_overview <- method_key %>%
    dplyr::arrange(.data$method, .data$model_name, .data$first_row) %>%
    dplyr::select(
      .data$method,
      .data$method_detail,
      .data$family,
      .data$model_name,
      .data$model,
      .data$residualizer,
      .data$residualizer_tag,
      .data$sem_c_order,
      .data$residualizer_c_order,
      .data$sem_exclusion,
      .data$residualizer_exclusion,
      .data$free_loadings,
      .data$bootstrap_B,
      .data$first_row,
      .data$uses_bootstrap_se,
      .data$has_sem_adjustment,
      .data$has_residualizer_adjustment,
      .data$has_confounder_adjustment,
      .data$n_replications,
      .data$N_values,
      .data$scenario_values
    ) %>%
    tibble::as_tibble()

  if (print_messages) {
    if (dplyr::n_distinct(df$scenario_id, na.rm = TRUE) > 1) {
      message("Multiple scenario_id values are present. Summaries keep scenarios separate.")
    }

    if (dplyr::n_distinct(df$N, na.rm = TRUE) > 1) {
      message("Multiple N values are present. Plots will separate sample sizes with facets.")
    }
  }

  list(
    df = df,
    method_key = method_key,
    method_levels = method_levels,
    n_levels = n_levels,
    palette = pal_method,
    model_overview = model_overview
  )
}

make_analysis_long <- function(
    df,
    effects = c("ARX", "ARY", "CXY", "CYX"),
    effect_map = default_effect_map,
    true_col_map = default_true_col_map
) {
  if (is.null(names(effect_map))) {
    stop("effect_map must be a named object, for example c(ARX = 'ARX', CXY = 'CXY').")
  }

  missing_effects <- setdiff(effects, names(effect_map))

  if (length(missing_effects) > 0) {
    stop(
      "These requested effects are not in effect_map: ",
      paste(missing_effects, collapse = ", ")
    )
  }

  get_one_mapping <- function(effect_name, map, map_name) {
    value <- map[[effect_name]]

    value <- unlist(value, recursive = TRUE, use.names = FALSE)
    value <- as.character(value)
    value <- value[!is.na(value) & value != ""]

    if (length(value) == 0) {
      stop(
        "Effect ", effect_name,
        " has no usable entry in ", map_name, "."
      )
    }

    if (length(value) > 1) {
      warning(
        "Effect ", effect_name,
        " maps to multiple columns in ", map_name,
        ": ",
        paste(value, collapse = ", "),
        ". Using the first one: ",
        value[1],
        call. = FALSE
      )
    }

    value[1]
  }

  get_column_or_na <- function(data, column_name) {
    column_name <- as.character(column_name)

    if (length(column_name) != 1L || is.na(column_name) || column_name == "") {
      return(rep(NA_real_, nrow(data)))
    }

    column_index <- match(column_name, names(data))

    if (is.na(column_index)) {
      return(rep(NA_real_, nrow(data)))
    }

    data[[column_index]]
  }

  pieces <- lapply(effects, function(effect_name) {
    estimate_col <- get_one_mapping(
      effect_name = effect_name,
      map = effect_map,
      map_name = "effect_map"
    )

    estimate_index <- match(estimate_col, names(df))

    if (is.na(estimate_index)) {
      return(NULL)
    }

    se_col <- paste0("se_", estimate_col)

    true_col <- if (estimate_col %in% names(true_col_map)) {
      get_one_mapping(
        effect_name = estimate_col,
        map = true_col_map,
        map_name = "true_col_map"
      )
    } else {
      NA_character_
    }

    estimate_value <- df[[estimate_index]]
    se_value <- get_column_or_na(df, se_col)
    true_value <- get_column_or_na(df, true_col)

    tibble::tibble(
      scenario_id = df$scenario_id,
      scenario_label = df$scenario_label,
      N = df$N,
      N_label = df$N_label,
      R = df$R,
      T = df$T,
      analysis_flag = df$analysis_flag,
      flag0 = df$flag0,
      flag1 = df$flag1,
      flag2 = df$flag2,
      flag3 = df$flag3,
      model_name = df$model_name,
      method_id = df$method_id,
      method = df$method,
      method_detail = df$method_detail,
      family = df$family,
      estimand = effect_name,
      estimate_col = estimate_col,
      estimate = estimate_value,
      se_estimate = se_value,
      true = true_value
    )
  })

  pieces <- pieces[!vapply(pieces, is.null, logical(1))]

  if (length(pieces) == 0) {
    stop("None of the requested estimate columns were found in results_df.")
  }

  dplyr::bind_rows(pieces) %>%
    dplyr::mutate(
      estimand = factor(.data$estimand, levels = effects),
      error = .data$estimate - .data$true,
      abs_error = abs(.data$error),
      sq_error = .data$error^2
    )
}


# ---- plotting helpers ----------------------------------------------------------------------------

add_standard_facets <- function(p, data, facet_row = NULL, scales = "fixed") {
  facet_cols <- character(0)

  if ("scenario_label" %in% names(data) && dplyr::n_distinct(data$scenario_label, na.rm = TRUE) > 1) {
    facet_cols <- c(facet_cols, "scenario_label")
  }

  if ("N_label" %in% names(data) && dplyr::n_distinct(data$N_label, na.rm = TRUE) > 1) {
    facet_cols <- c(facet_cols, "N_label")
  }

  row_part <- "."

  if (!is.null(facet_row) && facet_row %in% names(data) && dplyr::n_distinct(data[[facet_row]], na.rm = TRUE) > 1) {
    row_part <- facet_row
  }

  if (row_part == "." && length(facet_cols) == 0) {
    return(p)
  }

  col_part <- if (length(facet_cols) > 0) {
    paste(facet_cols, collapse = " + ")
  } else {
    "."
  }

  p + ggplot2::facet_grid(stats::as.formula(paste(row_part, "~", col_part)), scales = scales)
}

make_line_plot <- function(
    data,
    y,
    color_var,
    group_var,
    facet_row = NULL,
    y_se = NULL,
    palette = NULL,
    y_label = NULL,
    ref_line = NULL,
    ref_df = NULL,
    ref_y = NULL,
    ref_linetype = "dotted",
    zero_line = FALSE,
    percent_y = FALSE,
    clamp_01 = FALSE,
    scales = "fixed"
) {
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  data <- data %>%
    dplyr::filter(!is.na(.data[[y]]))

  if (nrow(data) == 0) {
    return(NULL)
  }

  p <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x = .data$T,
      y = .data[[y]],
      color = .data[[color_var]],
      group = .data[[group_var]]
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2)

  if (!is.null(y_se) && y_se %in% names(data)) {
    if (clamp_01) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = pmax(0, .data[[y]] - .data[[y_se]]),
            ymax = pmin(1, .data[[y]] + .data[[y_se]])
          ),
          width = 0.12,
          linewidth = 0.4
        )
    } else {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = .data[[y]] - .data[[y_se]],
            ymax = .data[[y]] + .data[[y_se]]
          ),
          width = 0.12,
          linewidth = 0.4
        )
    }
  }

  if (!is.null(ref_df) && !is.null(ref_y) && ref_y %in% names(ref_df) && nrow(ref_df) > 0) {
    p <- p +
      ggplot2::geom_line(
        data = ref_df,
        mapping = ggplot2::aes(x = .data$T, y = .data[[ref_y]]),
        inherit.aes = FALSE,
        color = "black",
        linetype = ref_linetype,
        linewidth = 0.7
      )
  }

  if (!is.null(ref_line)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
  }

  if (zero_line) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      )
  }

  p <- add_standard_facets(p, data, facet_row = facet_row, scales = scales)

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette, drop = FALSE)
  }

  p <- p + ggplot2::scale_x_continuous(breaks = sort(unique(data$T)))

  if (percent_y) {
    if (clamp_01) {
      p <- p +
        ggplot2::scale_y_continuous(
          labels = function(x) scales::percent(x, accuracy = 1),
          limits = c(0, 1)
        )
    } else {
      p <- p +
        ggplot2::scale_y_continuous(
          labels = function(x) scales::percent(x, accuracy = 1)
        )
    }
  }

  p +
    ggplot2::labs(
      x = "Occasion",
      y = y_label,
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 11),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10)
    )
}

make_family_plot <- function(
    data,
    y,
    family_name,
    method_key,
    palette,
    y_se = NULL,
    y_label = NULL,
    ref_line = NULL,
    zero_line = FALSE,
    percent_y = FALSE,
    clamp_01 = FALSE
) {
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  family_methods <- method_key %>%
    dplyr::filter(.data$family == .env$family_name) %>%
    dplyr::pull(.data$method) %>%
    as.character() %>%
    unique()

  family_data <- data %>%
    dplyr::filter(as.character(.data$method) %in% family_methods) %>%
    dplyr::mutate(method = factor(as.character(.data$method), levels = family_methods))

  if (nrow(family_data) == 0) {
    return(NULL)
  }

  family_palette <- palette[family_methods]

  make_line_plot(
    data = family_data,
    y = y,
    color_var = "method",
    group_var = "method_id",
    facet_row = "estimand",
    y_se = y_se,
    palette = family_palette,
    y_label = y_label,
    ref_line = ref_line,
    zero_line = zero_line,
    percent_y = percent_y,
    clamp_01 = clamp_01
  )
}

make_flag_bar_plot <- function(data, x, x_label, palette = NULL) {
  if (is.null(data) || nrow(data) == 0 || !x %in% names(data)) {
    return(NULL)
  }

  plot_df <- data %>%
    dplyr::filter(!is.na(.data[[x]])) %>%
    dplyr::arrange(dplyr::desc(.data[[x]])) %>%
    dplyr::mutate(method = factor(as.character(.data$method), levels = rev(unique(as.character(.data$method)))))

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[x]],
      y = .data$method,
      fill = .data$method
    )
  ) +
    ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
    ggplot2::scale_x_continuous(labels = function(v) scales::percent(v, accuracy = 1)) +
    ggplot2::labs(
      x = x_label,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 10)
    )

  p <- add_standard_facets(p, plot_df, facet_row = "family", scales = "free_y")

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette, drop = FALSE)
  }

  p
}

extract_true_delta_df <- function(data) {
  if (!"true_delta_t_vector" %in% names(data)) {
    return(tibble::tibble())
  }

  delta_base <- data %>%
    dplyr::select(
      dplyr::any_of(c("scenario_id", "scenario_label", "N", "N_label", "T")),
      .data$true_delta_t_vector
    ) %>%
    dplyr::filter(!vapply(.data$true_delta_t_vector, is.null, logical(1))) %>%
    dplyr::distinct()

  if (nrow(delta_base) == 0) {
    return(tibble::tibble())
  }

  delta_base %>%
    tidyr::unnest_longer(
      .data$true_delta_t_vector,
      values_to = "delta",
      indices_to = "delta_name"
    ) %>%
    tidyr::separate(
      .data$delta_name,
      into = c("outcome", "confounder"),
      sep = "__",
      remove = FALSE,
      fill = "right"
    ) %>%
    dplyr::mutate(
      outcome = sim_clean_chr(.data$outcome),
      confounder = sim_clean_chr(.data$confounder)
    ) %>%
    dplyr::arrange(.data$scenario_id, .data$N, .data$outcome, .data$confounder, .data$T)
}


# ---- main function -------------------------------------------------------------------------------

plot_overview_suite <- function(
    results_df,
    scenario_id = NULL,
    N = NULL,
    occasions = NULL,
    effects = c("ARX", "ARY", "CXY", "CYX"),
    effect_map = default_effect_map,
    true_col_map = default_true_col_map,
    keep_model_names = NULL,
    keep_families = NULL,
    method_order = NULL,
    n_order = NULL,
    n_labels = NULL,
    palette = NULL,
    drop_flagged = TRUE,
    included_analysis_flags = default_included_analysis_flags,
    alpha = 0.05,
    improper_top_n = 5,
    print_messages = TRUE
) {
  stopifnot(is.data.frame(results_df))

  if (!is.null(occasions) && !is.numeric(occasions)) {
    stop("occasions must be NULL or a numeric vector.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).")
  }

  if (!is.numeric(improper_top_n) || length(improper_top_n) != 1 || improper_top_n < 1) {
    stop("improper_top_n must be a single positive number.")
  }

  prep <- prepare_suite_data(
    results_df = results_df,
    scenario_id = scenario_id,
    N = N,
    occasions = occasions,
    keep_model_names = keep_model_names,
    keep_families = keep_families,
    method_order = method_order,
    n_order = n_order,
    n_labels = n_labels,
    palette = palette,
    print_messages = print_messages
  )

  df_all <- prep$df
  method_key <- prep$method_key
  method_levels <- prep$method_levels
  pal_method <- prep$palette

  analysis_long_all <- make_analysis_long(
    df = df_all,
    effects = effects,
    effect_map = effect_map,
    true_col_map = true_col_map
  )

  analysis_long_used <- analysis_long_all

  if (drop_flagged) {
    analysis_long_used <- analysis_long_used %>%
      dplyr::filter(.data$analysis_flag %in% included_analysis_flags)
  }

  summary_group_cols <- c(
    "scenario_id",
    "scenario_label",
    "N",
    "N_label",
    "T",
    "estimand",
    "method_id",
    "method",
    "method_detail",
    "family"
  )

  # ---- SEM performance summaries -----------------------------------------------------------------

  bias_df <- analysis_long_used %>%
    dplyr::filter(!is.na(.data$estimate), !is.na(.data$true)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(.data$true),
      mean_est = mean(.data$estimate, na.rm = TRUE),
      bias = .data$mean_est - .data$true,
      mcse_bias = stats::sd(.data$estimate, na.rm = TRUE) / sqrt(.data$nsim),
      rel_bias = dplyr::if_else(
        abs(.data$true) < 1e-12,
        NA_real_,
        (.data$mean_est - .data$true) / .data$true
      ),
      mcse_rel_bias = dplyr::if_else(
        abs(.data$true) < 1e-12,
        NA_real_,
        stats::sd(.data$estimate, na.rm = TRUE) / sqrt(.data$nsim) / abs(.data$true)
      ),
      .groups = "drop"
    )

  relbias_df <- bias_df

  rmse_df <- analysis_long_used %>%
    dplyr::filter(!is.na(.data$estimate), !is.na(.data$true)) %>%
    dplyr::mutate(sq_err = (.data$estimate - .data$true)^2) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_sq_err = mean(.data$sq_err, na.rm = TRUE),
      var_sq_err = stats::var(.data$sq_err, na.rm = TRUE),
      rmse = sqrt(.data$mean_sq_err),
      mcse_rmse = dplyr::if_else(
        .data$nsim <= 1 | is.na(.data$rmse) | .data$rmse == 0,
        NA_real_,
        sqrt(.data$var_sq_err / .data$nsim) / (2 * .data$rmse)
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(-.data$mean_sq_err, -.data$var_sq_err)

  se_df <- summarise_mean_mcse(
    data = analysis_long_used %>% dplyr::filter(!is.na(.data$se_estimate)),
    value_col = "se_estimate",
    mean_name = "mean_se",
    mcse_name = "mcse_mean_se",
    group_cols = summary_group_cols
  )

  crit_z <- stats::qnorm(1 - alpha / 2)

  detect_df <- analysis_long_used %>%
    dplyr::filter(!is.na(.data$estimate), !is.na(.data$se_estimate), .data$se_estimate > 0) %>%
    dplyr::mutate(
      z_value = .data$estimate / .data$se_estimate,
      detected = abs(.data$z_value) > crit_z
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(.data$true),
      detect_prob = mean(.data$detected, na.rm = TRUE),
      mcse_detect = sqrt(.data$detect_prob * (1 - .data$detect_prob) / .data$nsim),
      error_type = dplyr::if_else(abs(.data$true) < 1e-12, "Type I error", "Power"),
      .groups = "drop"
    )

  se_check <- analysis_long_used %>%
    dplyr::filter(!is.na(.data$estimate), !is.na(.data$se_estimate)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mc_sd = stats::sd(.data$estimate, na.rm = TRUE),
      mean_se = mean(.data$se_estimate, na.rm = TRUE),
      diff = .data$mean_se - .data$mc_sd,
      ratio = dplyr::if_else(.data$mc_sd > 0, .data$mean_se / .data$mc_sd, NA_real_),
      .groups = "drop"
    )

  # ---- ML summaries ------------------------------------------------------------------------------

  ml_metric_cols <- intersect(c("mse_x", "mse_y", "r2_x", "r2_y"), names(df_all))

  if (length(ml_metric_cols) > 0) {
    ml_base <- df_all

    if (drop_flagged) {
      ml_base <- ml_base %>%
        dplyr::filter(.data$analysis_flag %in% included_analysis_flags)
    }

    ml_long <- ml_base %>%
      dplyr::filter(.data$residualizer != "N", !is.na(.data$residualizer_tag), .data$residualizer_tag != "None") %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(ml_metric_cols),
        names_to = "metric_name",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        metric = dplyr::case_when(
          .data$metric_name %in% c("mse_x", "mse_y") ~ "MSE",
          .data$metric_name %in% c("r2_x", "r2_y") ~ "R2",
          TRUE ~ NA_character_
        ),
        target = dplyr::case_when(
          .data$metric_name %in% c("mse_x", "r2_x") ~ "X",
          .data$metric_name %in% c("mse_y", "r2_y") ~ "Y",
          TRUE ~ NA_character_
        ),
        ml_method = .data$residualizer_tag,
        ml_method_id = paste0("residualizer=", .data$residualizer)
      ) %>%
      dplyr::filter(!is.na(.data$metric), !is.na(.data$target), !is.na(.data$value))
  } else {
    ml_long <- tibble::tibble()
  }

  ml_group_cols <- c("scenario_id", "scenario_label", "N", "N_label", "T", "target", "ml_method_id", "ml_method")

  mse_df <- if (nrow(ml_long) > 0) {
    summarise_mean_mcse(
      data = ml_long %>% dplyr::filter(.data$metric == "MSE"),
      value_col = "value",
      mean_name = "mean_mse",
      mcse_name = "mcse_mean_mse",
      group_cols = ml_group_cols
    )
  } else {
    tibble::tibble()
  }

  r2_df <- if (nrow(ml_long) > 0) {
    summarise_mean_mcse(
      data = ml_long %>% dplyr::filter(.data$metric == "R2"),
      value_col = "value",
      mean_name = "mean_r2",
      mcse_name = "mcse_mean_r2",
      group_cols = ml_group_cols
    )
  } else {
    tibble::tibble()
  }

  theory_df <- if (all(c("true_r2_x", "true_r2_y") %in% names(df_all))) {
    df_all %>%
      dplyr::group_by(.data$scenario_id, .data$scenario_label, .data$N, .data$N_label, .data$T) %>%
      dplyr::summarise(
        true_r2_x = mean(.data$true_r2_x, na.rm = TRUE),
        true_r2_y = mean(.data$true_r2_y, na.rm = TRUE),
        theoretical_min_mse_x = 1 - .data$true_r2_x,
        theoretical_min_mse_y = 1 - .data$true_r2_y,
        .groups = "drop"
      )
  } else {
    tibble::tibble()
  }

  ml_levels <- infer_ordered_levels(
    c(as.character(mse_df$ml_method), as.character(r2_df$ml_method)),
    preferred = c("LM", "EN", "BCA", "XGB")
  )

  pal_ml <- make_palette(ml_levels)

  if (nrow(mse_df) > 0) {
    mse_df$ml_method <- factor(as.character(mse_df$ml_method), levels = ml_levels)
  }

  if (nrow(r2_df) > 0) {
    r2_df$ml_method <- factor(as.character(r2_df$ml_method), levels = ml_levels)
  }

  # ---- run-level summaries ------------------------------------------------------------------------

  run_level <- df_all %>%
    dplyr::group_by(
      .data$scenario_id,
      .data$scenario_label,
      .data$N,
      .data$N_label,
      .data$method_id,
      .data$method,
      .data$method_detail,
      .data$family,
      .data$R
    ) %>%
    dplyr::summarise(
      flag0 = dplyr::first(.data$flag0),
      flag1 = dplyr::first(.data$flag1),
      flag2 = dplyr::first(.data$flag2),
      flag3 = dplyr::first(.data$flag3),
      improper_reason = dplyr::first(.data$improper_reason),
      .groups = "drop"
    )

  flag_group_cols <- c("scenario_id", "scenario_label", "N", "N_label", "method_id", "method", "method_detail", "family")

  flag0_df <- run_level %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(flag_group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), prop_flag0 = mean(.data$flag0, na.rm = TRUE), .groups = "drop")

  flag1_df <- run_level %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(flag_group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), prop_flag1 = mean(.data$flag1, na.rm = TRUE), .groups = "drop")

  flag2_df <- run_level %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(flag_group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), prop_flag2 = mean(.data$flag2, na.rm = TRUE), .groups = "drop")

  flag3_df <- run_level %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(flag_group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), prop_flag3 = mean(.data$flag3, na.rm = TRUE), .groups = "drop")

  improper_reason_df <- run_level %>%
    dplyr::filter(.data$flag2 == 1 | .data$flag3 == 1 | !is.na(.data$improper_reason)) %>%
    dplyr::mutate(
      improper_reason = dplyr::if_else(
        is.na(.data$improper_reason),
        "unspecified improper reason",
        .data$improper_reason
      )
    ) %>%
    dplyr::group_by(
      .data$scenario_id,
      .data$scenario_label,
      .data$N,
      .data$N_label,
      .data$method_id,
      .data$method,
      .data$method_detail,
      .data$family,
      .data$improper_reason
    ) %>%
    dplyr::summarise(n_reason = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(
      n_improper_total = sum(.data$n_reason),
      prop_reason = .data$n_reason / .data$n_improper_total
    ) %>%
    dplyr::ungroup()

  improper_reason_top_df <- improper_reason_df %>%
    dplyr::group_by(.data$scenario_id, .data$N, .data$method_id, .data$method, .data$family) %>%
    dplyr::slice_max(order_by = .data$n_reason, n = improper_top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  true_delta_df <- extract_true_delta_df(df_all)

  # ---- factor ordering ----------------------------------------------------------------------------

  for (nm in c("method")) {
    relbias_df[[nm]] <- factor(as.character(relbias_df[[nm]]), levels = method_levels)
    rmse_df[[nm]] <- factor(as.character(rmse_df[[nm]]), levels = method_levels)
    se_df[[nm]] <- factor(as.character(se_df[[nm]]), levels = method_levels)
    detect_df[[nm]] <- factor(as.character(detect_df[[nm]]), levels = method_levels)
    se_check[[nm]] <- factor(as.character(se_check[[nm]]), levels = method_levels)
    flag0_df[[nm]] <- factor(as.character(flag0_df[[nm]]), levels = method_levels)
    flag1_df[[nm]] <- factor(as.character(flag1_df[[nm]]), levels = method_levels)
    flag2_df[[nm]] <- factor(as.character(flag2_df[[nm]]), levels = method_levels)
    flag3_df[[nm]] <- factor(as.character(flag3_df[[nm]]), levels = method_levels)
  }

  if (nrow(improper_reason_df) > 0) {
    improper_reason_df$method <- factor(as.character(improper_reason_df$method), levels = method_levels)
  }

  if (nrow(improper_reason_top_df) > 0) {
    improper_reason_top_df$method <- factor(as.character(improper_reason_top_df$method), levels = method_levels)
  }

  # ---- plots --------------------------------------------------------------------------------------

  plot_relbias_clpm <- make_family_plot(relbias_df, "rel_bias", "CLPM", method_key, pal_method, "mcse_rel_bias", "Relative bias", zero_line = TRUE)
  plot_relbias_riclpm <- make_family_plot(relbias_df, "rel_bias", "RICLPM", method_key, pal_method, "mcse_rel_bias", "Relative bias", zero_line = TRUE)
  plot_relbias_dpm <- make_family_plot(relbias_df, "rel_bias", "DPM", method_key, pal_method, "mcse_rel_bias", "Relative bias", zero_line = TRUE)

  plot_bias_clpm <- make_family_plot(bias_df, "bias", "CLPM", method_key, pal_method, "mcse_bias", "Bias", zero_line = TRUE)
  plot_bias_riclpm <- make_family_plot(bias_df, "bias", "RICLPM", method_key, pal_method, "mcse_bias", "Bias", zero_line = TRUE)
  plot_bias_dpm <- make_family_plot(bias_df, "bias", "DPM", method_key, pal_method, "mcse_bias", "Bias", zero_line = TRUE)

  plot_rmse_clpm <- make_family_plot(rmse_df, "rmse", "CLPM", method_key, pal_method, "mcse_rmse", "RMSE")
  plot_rmse_riclpm <- make_family_plot(rmse_df, "rmse", "RICLPM", method_key, pal_method, "mcse_rmse", "RMSE")
  plot_rmse_dpm <- make_family_plot(rmse_df, "rmse", "DPM", method_key, pal_method, "mcse_rmse", "RMSE")

  plot_se_clpm <- make_family_plot(se_df, "mean_se", "CLPM", method_key, pal_method, "mcse_mean_se", "Mean estimated SE")
  plot_se_riclpm <- make_family_plot(se_df, "mean_se", "RICLPM", method_key, pal_method, "mcse_mean_se", "Mean estimated SE")
  plot_se_dpm <- make_family_plot(se_df, "mean_se", "DPM", method_key, pal_method, "mcse_mean_se", "Mean estimated SE")

  plot_power_clpm <- make_family_plot(detect_df, "detect_prob", "CLPM", method_key, pal_method, "mcse_detect", "Detection probability", ref_line = alpha, percent_y = TRUE, clamp_01 = TRUE)
  plot_power_riclpm <- make_family_plot(detect_df, "detect_prob", "RICLPM", method_key, pal_method, "mcse_detect", "Detection probability", ref_line = alpha, percent_y = TRUE, clamp_01 = TRUE)
  plot_power_dpm <- make_family_plot(detect_df, "detect_prob", "DPM", method_key, pal_method, "mcse_detect", "Detection probability", ref_line = alpha, percent_y = TRUE, clamp_01 = TRUE)

  plot_se_ratio_clpm <- make_family_plot(se_check, "ratio", "CLPM", method_key, pal_method, NULL, "Estimated SE / Monte Carlo SD", ref_line = 1)
  plot_se_ratio_riclpm <- make_family_plot(se_check, "ratio", "RICLPM", method_key, pal_method, NULL, "Estimated SE / Monte Carlo SD", ref_line = 1)
  plot_se_ratio_dpm <- make_family_plot(se_check, "ratio", "DPM", method_key, pal_method, NULL, "Estimated SE / Monte Carlo SD", ref_line = 1)

  plot_se_diff_clpm <- make_family_plot(se_check, "diff", "CLPM", method_key, pal_method, NULL, "Estimated SE - Monte Carlo SD", ref_line = 0)
  plot_se_diff_riclpm <- make_family_plot(se_check, "diff", "RICLPM", method_key, pal_method, NULL, "Estimated SE - Monte Carlo SD", ref_line = 0)
  plot_se_diff_dpm <- make_family_plot(se_check, "diff", "DPM", method_key, pal_method, NULL, "Estimated SE - Monte Carlo SD", ref_line = 0)

  plot_se_ratio <- make_line_plot(
    data = se_check,
    y = "ratio",
    color_var = "method",
    group_var = "method_id",
    facet_row = "estimand",
    palette = pal_method,
    y_label = "Estimated SE / Monte Carlo SD",
    ref_line = 1
  )

  plot_se_diff <- make_line_plot(
    data = se_check,
    y = "diff",
    color_var = "method",
    group_var = "method_id",
    facet_row = "estimand",
    palette = pal_method,
    y_label = "Estimated SE - Monte Carlo SD",
    ref_line = 0
  )

  plot_mse_x <- make_line_plot(
    data = mse_df %>% dplyr::filter(.data$target == "X"),
    y = "mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF MSE",
    ref_df = theory_df,
    ref_y = "theoretical_min_mse_x"
  )

  plot_mse_y <- make_line_plot(
    data = mse_df %>% dplyr::filter(.data$target == "Y"),
    y = "mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF MSE",
    ref_df = theory_df,
    ref_y = "theoretical_min_mse_y"
  )

  plot_r2_x <- make_line_plot(
    data = r2_df %>% dplyr::filter(.data$target == "X"),
    y = "mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF R²",
    ref_df = theory_df,
    ref_y = "true_r2_x"
  )

  plot_r2_y <- make_line_plot(
    data = r2_df %>% dplyr::filter(.data$target == "Y"),
    y = "mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF R²",
    ref_df = theory_df,
    ref_y = "true_r2_y"
  )

  plot_mse <- make_line_plot(
    data = mse_df,
    y = "mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    facet_row = "target",
    palette = pal_ml,
    y_label = "Mean OOF MSE"
  )

  plot_r2 <- make_line_plot(
    data = r2_df,
    y = "mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    facet_row = "target",
    palette = pal_ml,
    y_label = "Mean OOF R²"
  )

  plot_true_delta <- NULL

  if (nrow(true_delta_df) > 0) {
    plot_true_delta <- ggplot2::ggplot(
      true_delta_df,
      ggplot2::aes(
        x = .data$T,
        y = .data$delta,
        color = .data$confounder,
        group = .data$confounder
      )
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::scale_x_continuous(breaks = sort(unique(true_delta_df$T))) +
      ggplot2::labs(
        x = "Occasion",
        y = "True delta",
        color = "Confounder"
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(legend.position = "bottom")

    plot_true_delta <- add_standard_facets(plot_true_delta, true_delta_df, facet_row = "outcome")
  }

  plot_flag0 <- make_flag_bar_plot(flag0_df, "prop_flag0", "Share flag0", pal_method)
  plot_flag1 <- make_flag_bar_plot(flag1_df, "prop_flag1", "Share flag1", pal_method)
  plot_flag2 <- make_flag_bar_plot(flag2_df, "prop_flag2", "Share flag2", pal_method)
  plot_flag3 <- make_flag_bar_plot(flag3_df, "prop_flag3", "Share flag3", pal_method)

  plot_improper_reasons <- NULL

  if (nrow(improper_reason_top_df) > 0) {
    plot_df <- improper_reason_top_df %>%
      dplyr::mutate(
        improper_reason = stats::reorder(.data$improper_reason, .data$prop_reason)
      )

    plot_improper_reasons <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data$prop_reason,
        y = .data$improper_reason,
        fill = .data$method
      )
    ) +
      ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
      ggplot2::facet_wrap(stats::as.formula("~ method"), scales = "free_y") +
      ggplot2::scale_x_continuous(labels = function(v) scales::percent(v, accuracy = 1)) +
      ggplot2::scale_fill_manual(values = pal_method, drop = FALSE) +
      ggplot2::labs(
        x = "Share within improper runs",
        y = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 9)
      )
  }

  plots <- list(
    plot_relbias_clpm = plot_relbias_clpm,
    plot_relbias_riclpm = plot_relbias_riclpm,
    plot_relbias_dpm = plot_relbias_dpm,
    plot_bias_clpm = plot_bias_clpm,
    plot_bias_riclpm = plot_bias_riclpm,
    plot_bias_dpm = plot_bias_dpm,
    plot_rmse_clpm = plot_rmse_clpm,
    plot_rmse_riclpm = plot_rmse_riclpm,
    plot_rmse_dpm = plot_rmse_dpm,
    plot_se_clpm = plot_se_clpm,
    plot_se_riclpm = plot_se_riclpm,
    plot_se_dpm = plot_se_dpm,
    plot_power_clpm = plot_power_clpm,
    plot_power_riclpm = plot_power_riclpm,
    plot_power_dpm = plot_power_dpm,
    plot_se_ratio_clpm = plot_se_ratio_clpm,
    plot_se_ratio_riclpm = plot_se_ratio_riclpm,
    plot_se_ratio_dpm = plot_se_ratio_dpm,
    plot_se_diff_clpm = plot_se_diff_clpm,
    plot_se_diff_riclpm = plot_se_diff_riclpm,
    plot_se_diff_dpm = plot_se_diff_dpm,
    plot_se_ratio = plot_se_ratio,
    plot_se_diff = plot_se_diff,
    plot_mse_x = plot_mse_x,
    plot_mse_y = plot_mse_y,
    plot_r2_x = plot_r2_x,
    plot_r2_y = plot_r2_y,
    plot_mse = plot_mse,
    plot_r2 = plot_r2,
    plot_true_delta = plot_true_delta,
    plot_flag0 = plot_flag0,
    plot_flag1 = plot_flag1,
    plot_flag2 = plot_flag2,
    plot_flag3 = plot_flag3,
    plot_improper_reasons = plot_improper_reasons
  )

  list(
    config = list(
      scenario_id = scenario_id,
      N = N,
      occasions = occasions,
      effects = effects,
      drop_flagged = drop_flagged,
      included_analysis_flags = included_analysis_flags,
      alpha = alpha
    ),

    raw_df = df_all,
    analysis_long = analysis_long_all,
    analysis_long_used = analysis_long_used,

    model_overview = prep$model_overview,
    method_key = method_key,
    method_levels = method_levels,
    palette = pal_method,

    bias_df = bias_df,
    relbias_df = relbias_df,
    rmse_df = rmse_df,
    se_df = se_df,
    detect_df = detect_df,
    se_check = se_check,

    mse_df = mse_df,
    r2_df = r2_df,
    theory_df = theory_df,
    true_delta_df = true_delta_df,

    flag0_df = flag0_df,
    flag1_df = flag1_df,
    flag2_df = flag2_df,
    flag3_df = flag3_df,
    improper_reason_df = improper_reason_df,

    plots = plots,

    plot_relbias_clpm = plot_relbias_clpm,
    plot_relbias_riclpm = plot_relbias_riclpm,
    plot_relbias_dpm = plot_relbias_dpm,
    plot_bias_clpm = plot_bias_clpm,
    plot_bias_riclpm = plot_bias_riclpm,
    plot_bias_dpm = plot_bias_dpm,
    plot_rmse_clpm = plot_rmse_clpm,
    plot_rmse_riclpm = plot_rmse_riclpm,
    plot_rmse_dpm = plot_rmse_dpm,
    plot_se_clpm = plot_se_clpm,
    plot_se_riclpm = plot_se_riclpm,
    plot_se_dpm = plot_se_dpm,
    plot_power_clpm = plot_power_clpm,
    plot_power_riclpm = plot_power_riclpm,
    plot_power_dpm = plot_power_dpm,
    plot_se_ratio_clpm = plot_se_ratio_clpm,
    plot_se_ratio_riclpm = plot_se_ratio_riclpm,
    plot_se_ratio_dpm = plot_se_ratio_dpm,
    plot_se_diff_clpm = plot_se_diff_clpm,
    plot_se_diff_riclpm = plot_se_diff_riclpm,
    plot_se_diff_dpm = plot_se_diff_dpm,
    plot_se_ratio = plot_se_ratio,
    plot_se_diff = plot_se_diff,
    plot_mse_x = plot_mse_x,
    plot_mse_y = plot_mse_y,
    plot_r2_x = plot_r2_x,
    plot_r2_y = plot_r2_y,
    plot_mse = plot_mse,
    plot_r2 = plot_r2,
    plot_true_delta = plot_true_delta,
    plot_flag0 = plot_flag0,
    plot_flag1 = plot_flag1,
    plot_flag2 = plot_flag2,
    plot_flag3 = plot_flag3,
    plot_improper_reasons = plot_improper_reasons
  )
}