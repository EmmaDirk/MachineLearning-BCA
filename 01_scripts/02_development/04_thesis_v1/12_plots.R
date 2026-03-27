# ============================================================
# libraries + setup
# ============================================================
library(here)
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))

# source plotting function
source(here("01_scripts", "02_development", "04_thesis_v1", "11_plotting.R"))

# ============================================================
# load data
# ============================================================
all_model_300_dfs <- readRDS(
  here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_300.rds"
  )
)

all_model_5000_dfs <- readRDS(
  here(
    "03_output", "02_thesis", "01_tests", "01_data",
    "02_constant_nl_omit",
    "00_all_10_model_dataframes_constant_2c_linear_5000.rds"
  )
)

all_model_300_df <- dplyr::bind_rows(all_model_300_dfs)
all_model_5000_df <- dplyr::bind_rows(all_model_5000_dfs)

str(all_model_300_df)
str(all_model_5000_df)

# ============================================================
# original plots
# ============================================================
out_300 <- plot_engine_results(all_model_300_df)
out_5000 <- plot_engine_results(all_model_5000_df)

# 1) plots for n = 300
out_300$plot_relbias
out_300$plot_rmse
out_300$plot_se

out_300$plot_relbias_clpm
out_300$plot_rmse_clpm
out_300$plot_se_clpm

out_300$plot_relbias_riclpm
out_300$plot_rmse_riclpm
out_300$plot_se_riclpm

out_300$plot_relbias_dpm
out_300$plot_rmse_dpm
out_300$plot_se_dpm

# 2) plots for n = 5000
out_5000$plot_relbias
out_5000$plot_rmse
out_5000$plot_se

out_5000$plot_relbias_clpm
out_5000$plot_rmse_clpm
out_5000$plot_se_clpm

out_5000$plot_relbias_riclpm
out_5000$plot_rmse_riclpm
out_5000$plot_se_riclpm

out_5000$plot_relbias_dpm
out_5000$plot_rmse_dpm
out_5000$plot_se_dpm

# ============================================================
# helpers for method labels
# ============================================================
model_map <- c(
  "C" = "CLPM",
  "R" = "RI-CLPM",
  "D" = "DPM"
)

resid_map <- c(
  "N" = "None",
  "none" = "None",
  "L" = "LM",
  "linear" = "LM",
  "X" = "XGB",
  "xgb" = "XGB"
)

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

make_method_label <- function(model, residualizer, exclusion, c_order,
                              free_loadings = NA, bootstrap_B = NA) {

  model_full <- dplyr::recode(
    as.character(model),
    !!!model_map,
    .default = as.character(model)
  )

  resid_full <- dplyr::recode(
    as.character(residualizer),
    !!!resid_map,
    .default = as.character(residualizer)
  )

  exclusion <- clean_text(exclusion)

  c_order_num <- suppressWarnings(as.numeric(c_order))
  free_loadings_num <- suppressWarnings(as.numeric(free_loadings))

  if (!is.na(free_loadings_num) && free_loadings_num == 1) {
    model_full <- dplyr::case_when(
      model_full == "RI-CLPM" ~ "fRI-CLPM",
      model_full == "DPM" ~ "fDPM",
      TRUE ~ model_full
    )
  }

  base_label <- dplyr::case_when(
    is.na(resid_full) ~ model_full,
    resid_full == "None" ~ model_full,
    TRUE ~ paste("BCA", resid_full, model_full)
  )

  is_unadjusted <- is.na(resid_full) || resid_full == "None"

  detail_part <- dplyr::case_when(
    is_unadjusted && model_full %in% c("RI-CLPM", "fRI-CLPM", "DPM", "fDPM") ~ NA_character_,
    is_unadjusted && model_full == "CLPM" && !is.na(c_order_num) ~ as.character(c_order_num),
    is_unadjusted && model_full == "CLPM" && is.na(c_order_num) ~ "0",
    !is.na(c_order_num) && !is.na(exclusion) ~ paste0(c_order_num, ":", exclusion),
    !is.na(c_order_num) ~ as.character(c_order_num),
    !is.na(exclusion) ~ exclusion,
    TRUE ~ NA_character_
  )

  if (is.na(detail_part)) {
    base_label
  } else {
    paste0(base_label, " (", detail_part, ")")
  }
}

make_method_id <- function(model, residualizer, exclusion, c_order,
                           free_loadings = NA, bootstrap_B = NA) {
  paste(
    paste0("model=", clean_text(model)),
    paste0("resid=", clean_text(residualizer)),
    paste0("excl=", clean_text(exclusion)),
    paste0("corder=", clean_text(c_order)),
    paste0("free=", clean_text(free_loadings)),
    sep = " | "
  )
}

# ============================================================
# main diagnostic function
# ============================================================
make_sim_diagnostics <- function(results_df, occasions = 2:5) {

  df_all <- results_df %>%
    dplyr::mutate(
      exclusion = clean_text(exclusion),
      c_order = clean_text(c_order),
      free_loadings = clean_text(free_loadings),
      bootstrap_B = clean_text(bootstrap_B),
      method = purrr::pmap_chr(
        list(model, residualizer, exclusion, c_order, free_loadings, bootstrap_B),
        make_method_label
      ),
      method_id = purrr::pmap_chr(
        list(model, residualizer, exclusion, c_order, free_loadings, bootstrap_B),
        make_method_id
      )
    ) %>%
    dplyr::filter(T %in% occasions)

  # ----------------------------------------------------------
  # 1. flag summaries
  # ----------------------------------------------------------
  flag_summary_method <- df_all %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_flag0 = sum(flag == 0, na.rm = TRUE),
      n_flag1 = sum(flag == 1, na.rm = TRUE),
      n_flag_other = sum(!is.na(flag) & !flag %in% c(0, 1)),
      prop_flagged = mean(flag != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prop_flagged), method)

  flag_summary_method_T <- df_all %>%
    dplyr::group_by(method, T) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_flag0 = sum(flag == 0, na.rm = TRUE),
      n_flag1 = sum(flag == 1, na.rm = TRUE),
      n_flag_other = sum(!is.na(flag) & !flag %in% c(0, 1)),
      prop_flagged = mean(flag != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prop_flagged), method, T)

  flag_summary_wide <- df_all %>%
    dplyr::count(method, T, flag) %>%
    dplyr::mutate(flag = paste0("flag_", flag)) %>%
    tidyr::pivot_wider(
      names_from = flag,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::arrange(method, T)

  # ----------------------------------------------------------
  # 2. long estimate and SE data before flag filtering
  # ----------------------------------------------------------
  est_long_all <- df_all %>%
    tidyr::pivot_longer(
      cols = c(ARX, ARY, CXY, CYX),
      names_to = "estimand",
      values_to = "estimate"
    )

  se_long_all <- df_all %>%
    tidyr::pivot_longer(
      cols = c(se_ARX, se_ARY, se_CXY, se_CYX),
      names_to = "se_name",
      values_to = "se_estimate"
    ) %>%
    dplyr::mutate(
      estimand = dplyr::recode(
        se_name,
        se_ARX = "ARX",
        se_ARY = "ARY",
        se_CXY = "CXY",
        se_CYX = "CYX"
      )
    )

  sim_long_all <- est_long_all %>%
    dplyr::left_join(
      se_long_all %>%
        dplyr::select(R, T, method_id, estimand, se_estimate),
      by = c("R", "T", "method_id", "estimand")
    )

  # ----------------------------------------------------------
  # 3. replication-loss diagnostics
  # ----------------------------------------------------------
  usable_summary <- sim_long_all %>%
    dplyr::group_by(method, method_id, T, estimand) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_flag0 = sum(flag == 0, na.rm = TRUE),
      n_est_missing = sum(is.na(estimate)),
      n_se_missing = sum(is.na(se_estimate)),
      n_both_present = sum(!is.na(estimate) & !is.na(se_estimate)),
      n_used_after_flag0 = sum(flag == 0 & !is.na(estimate) & !is.na(se_estimate)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(estimand, T, method)

  # ----------------------------------------------------------
  # 4. SE checks on usable replications only
  # ----------------------------------------------------------
  sim_long_used <- sim_long_all %>%
    dplyr::filter(flag == 0, !is.na(estimate), !is.na(se_estimate))

  se_check <- sim_long_used %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_est = mean(estimate, na.rm = TRUE),
      emp_sd = sd(estimate, na.rm = TRUE),
      mean_se = mean(se_estimate, na.rm = TRUE),
      median_se = median(se_estimate, na.rm = TRUE),
      rmse_se = sqrt(mean((se_estimate - emp_sd)^2, na.rm = TRUE)),
      diff = mean_se - emp_sd,
      rel_diff = dplyr::if_else(emp_sd > 0, (mean_se - emp_sd) / emp_sd, NA_real_),
      ratio = dplyr::if_else(emp_sd > 0, mean_se / emp_sd, NA_real_),
      mcse_mean_se = stats::sd(se_estimate, na.rm = TRUE) / sqrt(nsim),
      mcse_emp_sd = dplyr::if_else(nsim > 1, emp_sd / sqrt(2 * (nsim - 1)), NA_real_),
      z_gap = diff / sqrt(mcse_mean_se^2 + mcse_emp_sd^2),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      usable_summary %>%
        dplyr::select(method, method_id, T, estimand, n_total, n_flag0, n_used_after_flag0),
      by = c("method", "method_id", "T", "estimand")
    ) %>%
    dplyr::arrange(estimand, T, method)

  # ----------------------------------------------------------
  # 5. optional standardized-error check
  # ----------------------------------------------------------
  truth_map <- tibble::tibble(
    estimand = c("ARX", "ARY", "CXY", "CYX"),
    true_col = c("beta_x", "beta_y", "gamma_xy", "gamma_yx")
  )

  z_check <- sim_long_used %>%
    dplyr::left_join(truth_map, by = "estimand") %>%
    dplyr::mutate(
      true = dplyr::case_when(
        true_col == "beta_x" ~ beta_x,
        true_col == "beta_y" ~ beta_y,
        true_col == "gamma_xy" ~ gamma_xy,
        true_col == "gamma_yx" ~ gamma_yx,
        TRUE ~ NA_real_
      ),
      z = (estimate - true) / se_estimate
    ) %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_z = mean(z, na.rm = TRUE),
      sd_z = sd(z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(estimand, T, method)

  # ----------------------------------------------------------
  # 6. optional plots for SE calibration
  # ----------------------------------------------------------
  plot_se_ratio <- ggplot2::ggplot(
    se_check,
    ggplot2::aes(x = T, y = ratio, color = method, group = method_id)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ estimand, nrow = 1) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      x = "Occasion",
      y = "Mean reported SE / empirical SD",
      color = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10)
    )

  plot_se_diff <- ggplot2::ggplot(
    se_check,
    ggplot2::aes(x = T, y = diff, color = method, group = method_id)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ estimand, nrow = 1) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      x = "Occasion",
      y = "Mean reported SE - empirical SD",
      color = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10)
    )

  list(
    flag_summary_method = flag_summary_method,
    flag_summary_method_T = flag_summary_method_T,
    flag_summary_wide = flag_summary_wide,
    usable_summary = usable_summary,
    se_check = se_check,
    z_check = z_check,
    plot_se_ratio = plot_se_ratio,
    plot_se_diff = plot_se_diff
  )
}

# ============================================================
# run diagnostics for both sample sizes
# ============================================================
diag_300 <- make_sim_diagnostics(all_model_300_df, occasions = 2:5)
diag_5000 <- make_sim_diagnostics(all_model_5000_df, occasions = 2:5)

# ============================================================
# inspect results: n = 300
# ============================================================
diag_300$flag_summary_method
diag_300$flag_summary_method_T
diag_300$flag_summary_wide

diag_300$usable_summary
diag_300$se_check
diag_300$z_check

diag_300$plot_se_ratio
diag_300$plot_se_diff

# ============================================================
# inspect results: n = 5000
# ============================================================
diag_5000$flag_summary_method
diag_5000$flag_summary_method_T
diag_5000$flag_summary_wide

diag_5000$usable_summary
diag_5000$se_check
diag_5000$z_check

diag_5000$plot_se_ratio
diag_5000$plot_se_diff

# ============================================================
# optional: quickly inspect worst methods
# ============================================================
diag_300$se_check %>%
  dplyr::arrange(ratio) %>%
  print(n = 50)

diag_5000$se_check %>%
  dplyr::arrange(ratio) %>%
  print(n = 50)

diag_300$flag_summary_method_T %>%
  dplyr::arrange(dplyr::desc(prop_flagged)) %>%
  print(n = 50)

diag_5000$flag_summary_method_T %>%
  dplyr::arrange(dplyr::desc(prop_flagged)) %>%
  print(n = 50)