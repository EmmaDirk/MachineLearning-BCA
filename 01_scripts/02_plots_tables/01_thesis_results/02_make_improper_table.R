# =================================================================================================
#
# This script creates the final LaTeX table of improper-solution proportions from the combined
# simulation-results data frame.
#
# The input is the combined output produced by the simulation runner. It contains one row per
# replication, occasion, model, sample-size condition, and scenario.
#
# The table reports, for each selected scenario, sample size, and method, the fraction of replications
# whose main model fit was classified as mildly or severely improper.
# =================================================================================================

# ---- libraries ----------------------------------------------------------------------------------

library(here)
library(tidyverse)


# ---- project paths -------------------------------------------------------------------------------

# Define the project root.

study_root <- here::here()

# Define reusable directories.

data_dir <- file.path(
  study_root,
  "02_data",
  "01_thesis_results"
)

table_dir <- file.path(
  study_root,
  "03_output"
)

# Define input and output files.

combined_results_file <- file.path(
  data_dir,
  "s1234_N300_1000_2000.rds"
)

latex_file <- file.path(
  table_dir,
  "improper_fit_proportions.tex"
)


# ---- scenario and sample-size setup --------------------------------------------------------------

# The names define the scenario labels used in the LaTeX table.
# The values define the scenario_id values used inside the combined results data frame.

scenario_map <- c(
  `1` = 1L,
  `2` = 2L,
  `3` = 3L,
  `4` = 4L
)

n_order_for_table <- c(2000, 1000, 300)


# ---- table defaults ------------------------------------------------------------------------------

# These are the models shown in the final table.

default_keep_model_names <- c(
  "clpm_linear_confounders",
  "clpm_xgb_residualized",
  "riclpm_no_adjustment",
  "riclpm_linear_confounders",
  "riclpm_xgb_residualized"
)

# This order controls the method columns in the LaTeX table.

default_method_order <- c(
  "CLPM adj",
  "CLPM BCA",
  "RICLPM",
  "RICLPM adj",
  "RICLPM BCA"
)

keep_model_names_for_table <- default_keep_model_names
method_order_for_table <- default_method_order

latex_caption <- "Fraction of Improper Solutions by Sample Size, Method, and Scenario"
latex_label <- "tab:improper-all-scenarios"


# ---- load combined data --------------------------------------------------------------------------

dat_all <- readRDS(combined_results_file)


# ---- scenario lists ------------------------------------------------------------------------------

# Split the combined results data frame into one list per scenario.
# Each scenario list contains one data frame per requested sample size.

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


dat_by_scenario <- purrr::imap(
  scenario_map,
  ~ make_scenario_list_from_combined(
    dat_all = dat_all,
    scenario_value = .x,
    n_order = n_order_for_table
  )
)


# ---- small utility helpers -----------------------------------------------------------------------

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

safe_max_flag <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))

  if (all(is.na(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE)
  }
}

format_latex_number <- function(x) {
  sprintf("%.2f", x)
}


# ---- model and method labels ---------------------------------------------------------------------

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

      method = method_base
    ) %>%
    dplyr::select(-method_base)
}


# ---- fit-flag handling ---------------------------------------------------------------------------

# Improper analysis flags:
# - 2 = converged but mildly improper solution
# - 3 = converged but severely improper solution

improper_analysis_flags <- c(2L, 3L)


# ---- table data preparation ----------------------------------------------------------------------

prepare_table_data <- function(
    dat_list,
    keep_model_names = default_keep_model_names,
    method_order = default_method_order
) {

  stopifnot(is.list(dat_list))
  stopifnot(all(purrr::map_lgl(dat_list, is.data.frame)))

  df <- dplyr::bind_rows(dat_list)

  required_cols <- c(
    "R", "N", "analysis_flag",
    "model_name", "model", "residualizer",
    "sem_exclusion", "sem_c_order",
    "residualizer_exclusion", "residualizer_c_order",
    "free_loadings", "bootstrap_B"
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df <- df %>%
    dplyr::mutate(
      analysis_flag = int_or_na(analysis_flag),
      N = int_or_na(N),

      model_name = clean_chr(model_name),
      model = clean_chr(model),
      residualizer = clean_chr(residualizer),
      sem_exclusion = clean_chr(sem_exclusion),
      residualizer_exclusion = clean_chr(residualizer_exclusion),

      model_tag = model_label(model, free_loadings),
      residualizer_tag = resid_label(residualizer),

      sem_c_num = num_or_na(sem_c_order),
      resid_c_num = num_or_na(residualizer_c_order),

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

  df %>%
    dplyr::mutate(
      method = factor(
        method,
        levels = c(
          method_order[method_order %in% unique(method)],
          setdiff(unique(method), method_order)
        )
      )
    )
}


# ---- improper fit table --------------------------------------------------------------------------

make_improper_fit_table <- function(
    dat_list,
    n_order = c(2000, 1000, 300),
    keep_model_names = default_keep_model_names,
    method_order = default_method_order
) {

  prepare_table_data(
    dat_list = dat_list,
    keep_model_names = keep_model_names,
    method_order = method_order
  ) %>%
    dplyr::group_by(N, method, model_name, R) %>%
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


# ---- latex table builder -------------------------------------------------------------------------

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


# ---- create and export table ---------------------------------------------------------------------

improper_by_scenario <- purrr::map(
  dat_by_scenario,
  ~ make_improper_fit_table(
    dat_list = .x,
    n_order = n_order_for_table,
    keep_model_names = keep_model_names_for_table,
    method_order = method_order_for_table
  )
)

latex_table_all_scenarios <- make_latex_improper_table_all_scenarios(
  improper_list = improper_by_scenario,
  scenario_numbers = names(improper_by_scenario),
  n_order = n_order_for_table,
  method_cols = method_order_for_table,
  caption = latex_caption,
  label = latex_label
)

writeLines(
  latex_table_all_scenarios,
  con = latex_file
)