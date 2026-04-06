# This script takes the row-bound simulation output and builds publication-style plots
# for:
# - relative bias
# - RMSE
# - mean SE
#
# It also creates a clean method label from the compact simulation codes and returns
# the summarised plotting data for inspection.
# -------------------------------------------------------------------------------------------------

plot_engine_results <- function(results_df,
                                drop_flagged = TRUE,
                                occasions = 2:5) {

  # ------------------------------------------------------------
  # checks
  # ------------------------------------------------------------
  required_cols <- c(
    "R", "T", "analysis_flag", "flag0", "flag1", "flag2",
    "model", "residualizer", "exclusion", "c_order",
    "free_loadings", "bootstrap_B", "bootstrap_prop_success",
    "improper_reason", "bootstrap_issue_vector",
    "beta_x", "beta_y", "gamma_xy", "gamma_yx",
    "ARX", "ARY", "CXY", "CYX",
    "se_ARX", "se_ARY", "se_CXY", "se_CYX",
    "mse_x", "r2_x", "mse_y", "r2_y"
  )

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------
  model_map <- c(
    "C" = "CLPM",
    "R" = "RI-CLPM",
    "D" = "DPM"
  )

  resid_map <- c(
    "N" = "None",
    "L" = "LM",
    "X" = "XGB"
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

    # IMPORTANT:
    # keep the confounder order visible whenever adjustment is plotted,
    # keep CLPM (0) visible in the unadjusted case,
    # and omit brackets for simple unadjusted RI-CLPM / DPM variants
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

  truth_map <- tibble::tibble(
    estimand = c("ARX", "ARY", "CXY", "CYX"),
    true_col = c("beta_x", "beta_y", "gamma_xy", "gamma_yx")
  )

  se_map <- tibble::tibble(
    estimand = c("ARX", "ARY", "CXY", "CYX"),
    se_col = c("se_ARX", "se_ARY", "se_CXY", "se_CYX")
  )

  # ------------------------------------------------------------
  # prepare data
  # ------------------------------------------------------------
  df_all <- results_df %>%
    dplyr::mutate(
      exclusion = clean_text(exclusion),
      c_order = clean_text(c_order),
      free_loadings = clean_text(free_loadings),
      bootstrap_B = clean_text(bootstrap_B),
      improper_reason = clean_text(improper_reason),
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

  # classic performance plots optionally use only admissible main-analysis runs
  df <- df_all
  if (drop_flagged) {
    df <- df %>% dplyr::filter(analysis_flag == 0)
  }

  # ------------------------------------------------------------
  # long format for the classic performance plots
  # ------------------------------------------------------------
  long_est_df <- df %>%
    tidyr::pivot_longer(
      cols = c("ARX", "ARY", "CXY", "CYX"),
      names_to = "estimand",
      values_to = "estimate"
    ) %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::left_join(truth_map, by = "estimand") %>%
    dplyr::mutate(
      true = dplyr::case_when(
        true_col == "beta_x"   ~ beta_x,
        true_col == "beta_y"   ~ beta_y,
        true_col == "gamma_xy" ~ gamma_xy,
        true_col == "gamma_yx" ~ gamma_yx,
        TRUE ~ NA_real_
      )
    )

  long_se_df <- df %>%
    tidyr::pivot_longer(
      cols = c("se_ARX", "se_ARY", "se_CXY", "se_CYX"),
      names_to = "se_name",
      values_to = "se_estimate"
    ) %>%
    dplyr::filter(!is.na(se_estimate)) %>%
    dplyr::left_join(se_map, by = c("se_name" = "se_col")) %>%
    dplyr::filter(!is.na(estimand))

  # ------------------------------------------------------------
  # summaries for the classic performance plots
  # ------------------------------------------------------------
  relbias_df <- long_est_df %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      mean_est = mean(estimate, na.rm = TRUE),
      rel_bias = dplyr::if_else(
        is.na(true) | abs(true) < 1e-12,
        NA_real_,
        (mean_est - true) / true
      ),
      mcse_rel_bias = dplyr::if_else(
        is.na(true) | abs(true) < 1e-12,
        NA_real_,
        stats::sd(estimate, na.rm = TRUE) / sqrt(nsim) / abs(true)
      ),
      .groups = "drop"
    )

  rmse_df <- long_est_df %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      rmse = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
      mse_vals = list((estimate - true)^2),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mcse_rmse = {
        z <- unlist(mse_vals)
        n <- length(z)
        if (n <= 1 || is.na(rmse) || rmse == 0) {
          NA_real_
        } else {
          sqrt(stats::var(z) / n) / (2 * rmse)
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mse_vals)

  se_df <- long_se_df %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_se = mean(se_estimate, na.rm = TRUE),
      mcse_mean_se = stats::sd(se_estimate, na.rm = TRUE) / sqrt(nsim),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # SE diagnostics
  # these are computed on admissible main-analysis runs only
  # ------------------------------------------------------------
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

  sim_long_used <- sim_long_all %>%
    dplyr::filter(analysis_flag == 0, !is.na(estimate), !is.na(se_estimate))

  se_check <- sim_long_used %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      emp_sd = stats::sd(estimate, na.rm = TRUE),
      mean_se = mean(se_estimate, na.rm = TRUE),
      diff = mean_se - emp_sd,
      ratio = dplyr::if_else(emp_sd > 0, mean_se / emp_sd, NA_real_),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # flag summaries
  # ------------------------------------------------------------
  flag0_df <- df_all %>%
    dplyr::group_by(T, method_id, method) %>%
    dplyr::summarise(prop_flag = mean(flag0, na.rm = TRUE), .groups = "drop")

  flag1_df <- df_all %>%
    dplyr::group_by(T, method_id, method) %>%
    dplyr::summarise(prop_flag = mean(flag1, na.rm = TRUE), .groups = "drop")

  flag2_df <- df_all %>%
    dplyr::group_by(T, method_id, method) %>%
    dplyr::summarise(prop_flag = mean(flag2, na.rm = TRUE), .groups = "drop")

  # ------------------------------------------------------------
  # ML diagnostics
  # We keep X and Y separate, because the output stores them separately.
  # ------------------------------------------------------------
  ml_long_df <- df %>%
    tidyr::pivot_longer(
      cols = c(mse_x, r2_x, mse_y, r2_y),
      names_to = "ml_name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        ml_name %in% c("mse_x", "mse_y") ~ "MSE",
        ml_name %in% c("r2_x", "r2_y") ~ "R2",
        TRUE ~ NA_character_
      ),
      target = dplyr::case_when(
        ml_name %in% c("mse_x", "r2_x") ~ "X",
        ml_name %in% c("mse_y", "r2_y") ~ "Y",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(metric), !is.na(target), !is.na(value))

  mse_df <- ml_long_df %>%
    dplyr::filter(metric == "MSE") %>%
    dplyr::group_by(T, target, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_mse = mean(value, na.rm = TRUE),
      mcse_mean_mse = stats::sd(value, na.rm = TRUE) / sqrt(nsim),
      .groups = "drop"
    )

  r2_df <- ml_long_df %>%
    dplyr::filter(metric == "R2") %>%
    dplyr::group_by(T, target, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mean_r2 = mean(value, na.rm = TRUE),
      mcse_mean_r2 = stats::sd(value, na.rm = TRUE) / sqrt(nsim),
      .groups = "drop"
    )

  # average the X and Y stage-1 MSE within each replication so it can be linked
  # to one causal-bias quantity per replication and estimand
  ml_run_avg_df <- df %>%
    dplyr::transmute(
      R, T, method_id, method,
      mse_avg = rowMeans(cbind(mse_x, mse_y), na.rm = TRUE)
    )

  bias_mse_df <- long_est_df %>%
    dplyr::left_join(
      ml_run_avg_df,
      by = c("R", "T", "method_id", "method")
    ) %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      mean_mse = mean(mse_avg, na.rm = TRUE),
      rel_bias = {
        true_val <- dplyr::first(true)
        mean_est <- mean(estimate, na.rm = TRUE)
        if (is.na(true_val) || abs(true_val) < 1e-12) {
          NA_real_
        } else {
          (mean_est - true_val) / true_val
        }
      },
      .groups = "drop"
    )

  # keep X and Y separate for the R^2 vs MSE association plot
  r2_mse_df <- df %>%
    dplyr::select(R, T, method_id, method, mse_x, r2_x, mse_y, r2_y) %>%
    tidyr::pivot_longer(
      cols = c(mse_x, r2_x, mse_y, r2_y),
      names_to = "ml_name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        ml_name %in% c("mse_x", "mse_y") ~ "MSE",
        ml_name %in% c("r2_x", "r2_y") ~ "R2",
        TRUE ~ NA_character_
      ),
      target = dplyr::case_when(
        ml_name %in% c("mse_x", "r2_x") ~ "X",
        ml_name %in% c("mse_y", "r2_y") ~ "Y",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(metric), !is.na(target)) %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = value
    ) %>%
    dplyr::group_by(T, target, method_id, method) %>%
    dplyr::summarise(
      mean_mse = mean(MSE, na.rm = TRUE),
      mean_r2 = mean(R2, na.rm = TRUE),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # method ordering and palette
  # ------------------------------------------------------------
  method_key <- df_all %>%
    dplyr::distinct(method_id, method) %>%
    dplyr::mutate(
      family = dplyr::case_when(
        grepl("CLPM", method) ~ "CLPM",
        grepl("RI-CLPM", method) ~ "RI-CLPM",
        grepl("DPM", method) ~ "DPM",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::arrange(family, method)

  method_levels <- method_key$method

  pal_method <- setNames(
    scales::hue_pal()(length(method_levels)),
    method_levels
  )

  relbias_df$method <- factor(as.character(relbias_df$method), levels = method_levels)
  rmse_df$method <- factor(as.character(rmse_df$method), levels = method_levels)
  se_df$method <- factor(as.character(se_df$method), levels = method_levels)
  se_check$method <- factor(as.character(se_check$method), levels = method_levels)
  flag0_df$method <- factor(as.character(flag0_df$method), levels = method_levels)
  flag1_df$method <- factor(as.character(flag1_df$method), levels = method_levels)
  flag2_df$method <- factor(as.character(flag2_df$method), levels = method_levels)
  mse_df$method <- factor(as.character(mse_df$method), levels = method_levels)
  r2_df$method <- factor(as.character(r2_df$method), levels = method_levels)
  bias_mse_df$method <- factor(as.character(bias_mse_df$method), levels = method_levels)
  r2_mse_df$method <- factor(as.character(r2_mse_df$method), levels = method_levels)

  # ------------------------------------------------------------
  # plotting helper for the main metric plots with MCSE bars
  # ------------------------------------------------------------
  make_metric_plot <- function(data, metric, metric_se, add_zero_line = FALSE) {

    p <- data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = T,
          y = .data[[metric]],
          color = method,
          group = method_id
        )
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data[[metric]] - .data[[metric_se]],
          ymax = .data[[metric]] + .data[[metric_se]]
        ),
        width = 0.12,
        linewidth = 0.4
      ) +
      ggplot2::facet_wrap(~ estimand, nrow = 1, scales = "fixed") +
      ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$T))) +
      ggplot2::labs(
        title = NULL,
        x = "Occasion",
        y = NULL,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (add_zero_line) {
      p <- p + ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    p
  }

  # ------------------------------------------------------------
  # plotting helper for family-specific plots
  # ------------------------------------------------------------
  make_family_plot <- function(data, metric, metric_se, family_name,
                               add_zero_line = FALSE) {

    family_methods <- method_key %>%
      dplyr::filter(family == family_name) %>%
      dplyr::pull(method) %>%
      unique()

    family_data <- data %>%
      dplyr::filter(method %in% family_methods) %>%
      dplyr::mutate(method = factor(as.character(method), levels = family_methods))

    family_pal <- pal_method[family_methods]

    p <- family_data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = T,
          y = .data[[metric]],
          color = method,
          group = method_id
        )
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data[[metric]] - .data[[metric_se]],
          ymax = .data[[metric]] + .data[[metric_se]]
        ),
        width = 0.12,
        linewidth = 0.4
      ) +
      ggplot2::facet_wrap(~ estimand, nrow = 1, scales = "fixed") +
      ggplot2::scale_color_manual(values = family_pal, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = sort(unique(family_data$T))) +
      ggplot2::labs(
        title = NULL,
        x = "Occasion",
        y = NULL,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (add_zero_line) {
      p <- p + ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    p
  }

  # ------------------------------------------------------------
  # plotting helper for simple combined diagnostic plots
  # ------------------------------------------------------------
  make_simple_diagnostic_plot <- function(data, metric, facet_var, y_label,
                                          ref_line = NULL, percent_y = FALSE) {

    p <- data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = T,
          y = .data[[metric]],
          color = method,
          group = method_id
        )
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), nrow = 1) +
      ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$T))) +
      ggplot2::labs(
        title = NULL,
        x = "Occasion",
        y = y_label,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (!is.null(ref_line)) {
      p <- p + ggplot2::geom_hline(
        yintercept = ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    if (percent_y) {
      p <- p + ggplot2::scale_y_continuous(
        labels = function(x) scales::percent(x, accuracy = 1),
        limits = c(0, 1)
      )
    }

    p
  }

  # ------------------------------------------------------------
  # plotting helper for ML metric line plots with X/Y facets
  # ------------------------------------------------------------
  make_ml_line_plot <- function(data, metric, metric_se, y_label,
                                ref_line = NULL, percent_y = FALSE) {

    p <- data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = T,
          y = .data[[metric]],
          color = method,
          group = method_id
        )
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data[[metric]] - .data[[metric_se]],
          ymax = .data[[metric]] + .data[[metric_se]]
        ),
        width = 0.12,
        linewidth = 0.4
      ) +
      ggplot2::facet_wrap(~ target, nrow = 1, scales = "fixed") +
      ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$T))) +
      ggplot2::labs(
        title = NULL,
        x = "Occasion",
        y = y_label,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (!is.null(ref_line)) {
      p <- p + ggplot2::geom_hline(
        yintercept = ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    if (percent_y) {
      p <- p + ggplot2::scale_y_continuous(
        labels = function(x) scales::percent(x, accuracy = 1),
        limits = c(0, 1)
      )
    }

    p
  }

  # ------------------------------------------------------------
  # plotting helper for scatter-type association plots
  # ------------------------------------------------------------
  make_scatter_diagnostic_plot <- function(data, x_var, y_var, facet_var,
                                           x_label, y_label,
                                           x_ref_line = NULL, y_ref_line = NULL,
                                           percent_x = FALSE, percent_y = FALSE) {

    p <- data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data[[x_var]],
          y = .data[[y_var]],
          color = method,
          group = method_id
        )
      ) +
      ggplot2::geom_path(linewidth = 0.7, alpha = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(shape = factor(T)),
        size = 2.2
      ) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), nrow = 1, scales = "free") +
      ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
      ggplot2::labs(
        title = NULL,
        x = x_label,
        y = y_label,
        color = NULL,
        shape = "Occasion"
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (!is.null(x_ref_line)) {
      p <- p + ggplot2::geom_vline(
        xintercept = x_ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    if (!is.null(y_ref_line)) {
      p <- p + ggplot2::geom_hline(
        yintercept = y_ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    if (percent_x) {
      p <- p + ggplot2::scale_x_continuous(
        labels = function(x) scales::percent(x, accuracy = 1)
      )
    }

    if (percent_y) {
      p <- p + ggplot2::scale_y_continuous(
        labels = function(x) scales::percent(x, accuracy = 1)
      )
    }

    p
  }

  # ------------------------------------------------------------
  # classic plots
  # ------------------------------------------------------------
  plot_relbias <- make_metric_plot(
    data = relbias_df,
    metric = "rel_bias",
    metric_se = "mcse_rel_bias",
    add_zero_line = TRUE
  )

  plot_rmse <- make_metric_plot(
    data = rmse_df,
    metric = "rmse",
    metric_se = "mcse_rmse",
    add_zero_line = FALSE
  )

  plot_se <- make_metric_plot(
    data = se_df,
    metric = "mean_se",
    metric_se = "mcse_mean_se",
    add_zero_line = FALSE
  )

  plot_relbias_clpm <- make_family_plot(
    data = relbias_df,
    metric = "rel_bias",
    metric_se = "mcse_rel_bias",
    family_name = "CLPM",
    add_zero_line = TRUE
  )

  plot_relbias_riclpm <- make_family_plot(
    data = relbias_df,
    metric = "rel_bias",
    metric_se = "mcse_rel_bias",
    family_name = "RI-CLPM",
    add_zero_line = TRUE
  )

  plot_relbias_dpm <- make_family_plot(
    data = relbias_df,
    metric = "rel_bias",
    metric_se = "mcse_rel_bias",
    family_name = "DPM",
    add_zero_line = TRUE
  )

  plot_rmse_clpm <- make_family_plot(
    data = rmse_df,
    metric = "rmse",
    metric_se = "mcse_rmse",
    family_name = "CLPM",
    add_zero_line = FALSE
  )

  plot_rmse_riclpm <- make_family_plot(
    data = rmse_df,
    metric = "rmse",
    metric_se = "mcse_rmse",
    family_name = "RI-CLPM",
    add_zero_line = FALSE
  )

  plot_rmse_dpm <- make_family_plot(
    data = rmse_df,
    metric = "rmse",
    metric_se = "mcse_rmse",
    family_name = "DPM",
    add_zero_line = FALSE
  )

  plot_se_clpm <- make_family_plot(
    data = se_df,
    metric = "mean_se",
    metric_se = "mcse_mean_se",
    family_name = "CLPM",
    add_zero_line = FALSE
  )

  plot_se_riclpm <- make_family_plot(
    data = se_df,
    metric = "mean_se",
    metric_se = "mcse_mean_se",
    family_name = "RI-CLPM",
    add_zero_line = FALSE
  )

  plot_se_dpm <- make_family_plot(
    data = se_df,
    metric = "mean_se",
    metric_se = "mcse_mean_se",
    family_name = "DPM",
    add_zero_line = FALSE
  )

  # ------------------------------------------------------------
  # combined SE diagnostic plots
  # ------------------------------------------------------------
  plot_se_ratio <- make_simple_diagnostic_plot(
    data = se_check,
    metric = "ratio",
    facet_var = "estimand",
    y_label = "Mean reported SE / empirical SD",
    ref_line = 1,
    percent_y = FALSE
  )

  plot_se_diff <- make_simple_diagnostic_plot(
    data = se_check,
    metric = "diff",
    facet_var = "estimand",
    y_label = "Mean reported SE - empirical SD",
    ref_line = 0,
    percent_y = FALSE
  )

  # ------------------------------------------------------------
  # new ML plots
  # ------------------------------------------------------------
  plot_mse <- make_ml_line_plot(
    data = mse_df,
    metric = "mean_mse",
    metric_se = "mcse_mean_mse",
    y_label = "Mean OOF MSE",
    ref_line = NULL,
    percent_y = FALSE
  )

  plot_r2 <- make_ml_line_plot(
    data = r2_df,
    metric = "mean_r2",
    metric_se = "mcse_mean_r2",
    y_label = "Mean OOF R^2",
    ref_line = 0,
    percent_y = FALSE
  )

  plot_bias_vs_mse <- make_scatter_diagnostic_plot(
    data = bias_mse_df,
    x_var = "mean_mse",
    y_var = "rel_bias",
    facet_var = "estimand",
    x_label = "Mean OOF MSE",
    y_label = "Relative bias",
    x_ref_line = NULL,
    y_ref_line = 0,
    percent_x = FALSE,
    percent_y = FALSE
  )

  plot_r2_vs_mse <- make_scatter_diagnostic_plot(
    data = r2_mse_df,
    x_var = "mean_mse",
    y_var = "mean_r2",
    facet_var = "target",
    x_label = "Mean OOF MSE",
    y_label = "Mean OOF R^2",
    x_ref_line = NULL,
    y_ref_line = 0,
    percent_x = FALSE,
    percent_y = FALSE
  )

  # ------------------------------------------------------------
  # separate flag plots
  # ------------------------------------------------------------
  plot_flag0 <- ggplot2::ggplot(
    flag0_df,
    ggplot2::aes(x = T, y = prop_flag, color = method, group = method_id)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(flag0_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Proportion Flag 0",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 12)
    )

  plot_flag1 <- ggplot2::ggplot(
    flag1_df,
    ggplot2::aes(x = T, y = prop_flag, color = method, group = method_id)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(flag1_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Proportion Flag 1",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 12)
    )

  plot_flag2 <- ggplot2::ggplot(
    flag2_df,
    ggplot2::aes(x = T, y = prop_flag, color = method, group = method_id)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(flag2_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Proportion Flag 2",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 12)
    )

  # ------------------------------------------------------------
  # return everything for downstream use
  # ------------------------------------------------------------
  list(
    # summarised data
    relbias_df = relbias_df,
    rmse_df = rmse_df,
    se_df = se_df,
    se_check = se_check,
    flag0_df = flag0_df,
    flag1_df = flag1_df,
    flag2_df = flag2_df,
    mse_df = mse_df,
    r2_df = r2_df,
    bias_mse_df = bias_mse_df,
    r2_mse_df = r2_mse_df,

    # classic plots
    plot_relbias = plot_relbias,
    plot_rmse = plot_rmse,
    plot_se = plot_se,
    plot_relbias_clpm = plot_relbias_clpm,
    plot_relbias_riclpm = plot_relbias_riclpm,
    plot_relbias_dpm = plot_relbias_dpm,
    plot_rmse_clpm = plot_rmse_clpm,
    plot_rmse_riclpm = plot_rmse_riclpm,
    plot_rmse_dpm = plot_rmse_dpm,
    plot_se_clpm = plot_se_clpm,
    plot_se_riclpm = plot_se_riclpm,
    plot_se_dpm = plot_se_dpm,

    # SE diagnostic plots
    plot_se_ratio = plot_se_ratio,
    plot_se_diff = plot_se_diff,

    # new ML plots
    plot_mse = plot_mse,
    plot_r2 = plot_r2,
    plot_bias_vs_mse = plot_bias_vs_mse,
    plot_r2_vs_mse = plot_r2_vs_mse,

    # flag plots
    plot_flag0 = plot_flag0,
    plot_flag1 = plot_flag1,
    plot_flag2 = plot_flag2
  )
}