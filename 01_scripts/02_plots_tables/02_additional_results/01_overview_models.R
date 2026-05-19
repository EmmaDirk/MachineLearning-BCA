# =================================================================================================
#
# Helper functions for combined simulation-result data frames.
#
# Main function:
#   overview_models()
#
# Expected input:
#   One combined post-processing data frame with one row per replication, occasion, model,
#   scenario, and sample-size condition.
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
      dplyr::filter(.data$scenario_id %in% scenario_id)
  }

  if (!is.null(N)) {
    if (!"N" %in% names(df)) {
      stop("N was supplied, but results_df has no N column.")
    }

    df <- df %>%
      dplyr::filter(.data$N %in% N)
  }

  if (!is.null(occasions)) {
    if (!"T" %in% names(df)) {
      stop("occasions was supplied, but results_df has no T column.")
    }

    df <- df %>%
      dplyr::filter(.data$T %in% occasions)
  }

  if (!is.null(keep_model_names)) {
    if (!"model_name" %in% names(df)) {
      stop("keep_model_names was supplied, but results_df has no model_name column.")
    }

    df <- df %>%
      dplyr::filter(.data$model_name %in% keep_model_names)
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
      dplyr::filter(.data$model %in% wanted_codes)
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
    dplyr::select(-.data$method_base)

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
