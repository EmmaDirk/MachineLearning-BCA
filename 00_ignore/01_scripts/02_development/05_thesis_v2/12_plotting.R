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
    "R", "T", "flag0", "flag1", "flag2",
    "model", "residualizer", "exclusion", "c_order",
    "beta_x", "beta_y", "gamma_xy", "gamma_yx",
    "ARX", "ARY", "CXY", "CYX",
    "se_ARX", "se_ARY", "se_CXY", "se_CYX"
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
  if (!("free_loadings" %in% names(results_df))) {
    results_df$free_loadings <- NA_integer_
  }

  if (!("bootstrap_B" %in% names(results_df))) {
    results_df$bootstrap_B <- NA_integer_
  }

  if (!("bootstrap_prop_success" %in% names(results_df))) {
    results_df$bootstrap_prop_success <- NA_real_
  }

  if (!("improper_reason" %in% names(results_df))) {
    results_df$improper_reason <- NA_character_
  }

  if (!("bootstrap_issue_vector" %in% names(results_df))) {
    results_df$bootstrap_issue_vector <- rep(list(NA_character_), nrow(results_df))
  }

  # retain backward compatibility: if analysis_flag is absent, reconstruct it
  # from one-hot flag columns only when they are one-hot
  if (!("analysis_flag" %in% names(results_df))) {
    results_df$analysis_flag <- dplyr::case_when(
      results_df$flag0 == 1 ~ 0L,
      results_df$flag1 == 1 ~ 1L,
      results_df$flag2 == 1 ~ 2L,
      TRUE ~ NA_integer_
    )
  }

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
  # flag plots
  # IMPORTANT:
  # for one-step methods, flag0/1/2 are one-hot at the replication level
  # for bootstrap methods, flag0/1/2 are already bootstrap proportions
  # so taking the mean here gives the correct method-level percentages
  # ------------------------------------------------------------
  flag_plot_df <- df_all %>%
    dplyr::group_by(T, method_id, method) %>%
    dplyr::summarise(
      flag0 = mean(flag0, na.rm = TRUE),
      flag1 = mean(flag1, na.rm = TRUE),
      flag2 = mean(flag2, na.rm = TRUE),
      .groups = "drop"
    )

  flag0_df <- flag_plot_df %>%
    dplyr::transmute(T, method_id, method, prop_flag = flag0)

  flag1_df <- flag_plot_df %>%
    dplyr::transmute(T, method_id, method, prop_flag = flag1)

  flag2_df <- flag_plot_df %>%
    dplyr::transmute(T, method_id, method, prop_flag = flag2)

  # ------------------------------------------------------------
  # bootstrap success plot
  # ------------------------------------------------------------
  bootstrap_df <- df_all %>%
    dplyr::filter(!is.na(bootstrap_prop_success)) %>%
    dplyr::group_by(T, method_id, method) %>%
    dplyr::summarise(
      mean_bootstrap_prop_success = mean(bootstrap_prop_success, na.rm = TRUE),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # extra diagnostics: method-level summaries
  # IMPORTANT:
  # for bootstrap methods, mean_flag0/1/2 are mean bootstrap proportions
  # not proportions of main-analysis runs
  # prop_main_flagged still refers only to the main fit
  # ------------------------------------------------------------
  flag_summary_method <- df_all %>%
    dplyr::group_by(method, method_id) %>%
    dplyr::summarise(
      n_total = dplyr::n_distinct(R),
      mean_flag0 = mean(flag0, na.rm = TRUE),
      mean_flag1 = mean(flag1, na.rm = TRUE),
      mean_flag2 = mean(flag2, na.rm = TRUE),
      prop_main_flagged = mean(analysis_flag != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(prop_main_flagged), method)

  flag_summary_method_T <- df_all %>%
    dplyr::group_by(method, method_id, T) %>%
    dplyr::summarise(
      n_total = dplyr::n_distinct(R),
      mean_flag0 = mean(flag0, na.rm = TRUE),
      mean_flag1 = mean(flag1, na.rm = TRUE),
      mean_flag2 = mean(flag2, na.rm = TRUE),
      prop_main_flagged = mean(analysis_flag != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(method, T)

  # ------------------------------------------------------------
  # improper reason summaries for the main fit
  # ------------------------------------------------------------
  improper_reason_df <- df_all %>%
    dplyr::filter(analysis_flag == 2, !is.na(improper_reason)) %>%
    dplyr::group_by(method, method_id, T, improper_reason) %>%
    dplyr::summarise(
      n = dplyr::n_distinct(R),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(prop_within_method_T = n / sum(n)) %>%
    dplyr::ungroup()

  improper_reason_top_df <- improper_reason_df %>%
    dplyr::group_by(method, method_id, improper_reason) %>%
    dplyr::summarise(
      n = sum(n),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(prop_within_method = n / sum(n)) %>%
    dplyr::ungroup()

  improper_reason_plot_df <- improper_reason_df %>%
    dplyr::group_by(method, method_id, T) %>%
    dplyr::slice_max(order_by = prop_within_method_T, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # ------------------------------------------------------------
  # bootstrap issue summaries
  # ------------------------------------------------------------
  bootstrap_issue_long <- df_all %>%
    dplyr::filter(!is.na(bootstrap_prop_success)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      bootstrap_issue_vector_flat = list(
        if (is.list(bootstrap_issue_vector)) {
          unlist(bootstrap_issue_vector, use.names = FALSE)
        } else {
          as.character(bootstrap_issue_vector)
        }
      )
    ) %>%
    tidyr::unnest_longer(bootstrap_issue_vector_flat, values_to = "bootstrap_issue") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      bootstrap_issue = clean_text(bootstrap_issue)
    ) %>%
    dplyr::filter(!is.na(bootstrap_issue))

  bootstrap_issue_df <- bootstrap_issue_long %>%
    dplyr::group_by(method, method_id, T, bootstrap_issue) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(prop_within_method_T = n / sum(n)) %>%
    dplyr::ungroup()

  bootstrap_issue_top_df <- bootstrap_issue_df %>%
    dplyr::group_by(method, method_id, bootstrap_issue) %>%
    dplyr::summarise(
      n = sum(n),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(prop_within_method = n / sum(n)) %>%
    dplyr::ungroup()

  bootstrap_issue_plot_df <- bootstrap_issue_df %>%
    dplyr::filter(bootstrap_issue != "proper") %>%
    dplyr::group_by(method, method_id, T) %>%
    dplyr::slice_max(order_by = prop_within_method_T, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # ------------------------------------------------------------
  # ordering
  # ------------------------------------------------------------
  method_key <- df_all %>%
    dplyr::distinct(
      model, residualizer, exclusion, c_order,
      free_loadings, method_id, method
    ) %>%
    dplyr::mutate(
      family = dplyr::case_when(
        model == "C" ~ "CLPM",
        model == "R" ~ "RI-CLPM",
        model == "D" ~ "DPM",
        TRUE ~ "Other"
      ),
      family_order = dplyr::case_when(
        family == "CLPM" ~ 1L,
        family == "RI-CLPM" ~ 2L,
        family == "DPM" ~ 3L,
        TRUE ~ 99L
      )
    ) %>%
    dplyr::arrange(family_order, method)

  method_levels <- method_key %>%
    dplyr::pull(method) %>%
    unique()

  relbias_df <- relbias_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  rmse_df <- rmse_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  se_df <- se_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  se_check <- se_check %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  flag0_df <- flag0_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  flag1_df <- flag1_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  flag2_df <- flag2_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  bootstrap_df <- bootstrap_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  improper_reason_plot_df <- improper_reason_plot_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  bootstrap_issue_plot_df <- bootstrap_issue_plot_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  # ------------------------------------------------------------
  # palette
  # ------------------------------------------------------------
  make_family_palette <- function(methods, start_col, end_col) {
    n <- length(methods)

    if (n == 0) {
      return(stats::setNames(character(0), character(0)))
    }

    if (n == 1) {
      cols <- start_col
    } else {
      cols <- grDevices::colorRampPalette(c(start_col, end_col))(n)
    }

    stats::setNames(cols, methods)
  }

  clpm_methods <- method_key %>%
    dplyr::filter(family == "CLPM") %>%
    dplyr::pull(method)

  riclpm_methods <- method_key %>%
    dplyr::filter(family == "RI-CLPM") %>%
    dplyr::pull(method)

  dpm_methods <- method_key %>%
    dplyr::filter(family == "DPM") %>%
    dplyr::pull(method)

  pal_method <- c(
    make_family_palette(clpm_methods, "#39FF14", "#0033FF"),
    make_family_palette(riclpm_methods, "#FFD600", "#FF1F1F"),
    make_family_palette(dpm_methods, "#FF1493", "#7A00FF")
  )

  pal_method <- pal_method[method_levels]

  # ------------------------------------------------------------
  # plotting helper for the classic performance plots
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
  # family-specific plotting helper
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
      legend.text = ggplot2::element_text(size = 11)
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
      legend.text = ggplot2::element_text(size = 11)
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
      legend.text = ggplot2::element_text(size = 11)
    )

  # ------------------------------------------------------------
  # bootstrap success plot
  # ------------------------------------------------------------
  plot_bootstrap_prop_success <- ggplot2::ggplot(
    bootstrap_df,
    ggplot2::aes(
      x = T,
      y = mean_bootstrap_prop_success,
      color = method,
      group = method_id
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(bootstrap_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Bootstrap proportion success",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11)
    )

  # ------------------------------------------------------------
  # improper-reason plot
  # ------------------------------------------------------------
  plot_improper_reason_top <- ggplot2::ggplot(
    improper_reason_plot_df,
    ggplot2::aes(
      x = T,
      y = prop_within_method_T,
      color = method,
      group = method_id
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(
      ggplot2::aes(label = improper_reason),
      size = 3,
      hjust = 0,
      nudge_x = 0.05,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(improper_reason_plot_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Top improper-reason share",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11)
    )

  # ------------------------------------------------------------
  # bootstrap-issue plot
  # ------------------------------------------------------------
  plot_bootstrap_issue_top <- ggplot2::ggplot(
    bootstrap_issue_plot_df,
    ggplot2::aes(
      x = T,
      y = prop_within_method_T,
      color = method,
      group = method_id
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(
      ggplot2::aes(label = bootstrap_issue),
      size = 3,
      hjust = 0,
      nudge_x = 0.05,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal_method, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = sort(unique(bootstrap_issue_plot_df$T))) +
    ggplot2::scale_y_continuous(
      labels = function(x) scales::percent(x, accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Occasion",
      y = "Top bootstrap-issue share",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11)
    )

  # ------------------------------------------------------------
  # method key
  # ------------------------------------------------------------
  method_key <- method_key %>%
    dplyr::select(
      model, residualizer, exclusion, c_order,
      free_loadings, method_id, method, family
    )

  # ------------------------------------------------------------
  # return
  # ------------------------------------------------------------
  list(
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

    plot_se_ratio = plot_se_ratio,
    plot_se_diff = plot_se_diff,

    plot_flag0 = plot_flag0,
    plot_flag1 = plot_flag1,
    plot_flag2 = plot_flag2,

    plot_bootstrap_prop_success = plot_bootstrap_prop_success,
    plot_improper_reason_top = plot_improper_reason_top,
    plot_bootstrap_issue_top = plot_bootstrap_issue_top,

    relbias_df = relbias_df,
    rmse_df = rmse_df,
    se_df = se_df,
    se_check = se_check,

    flag_summary_method = flag_summary_method,
    flag_summary_method_T = flag_summary_method_T,

    improper_reason_df = improper_reason_df,
    improper_reason_top_df = improper_reason_top_df,

    bootstrap_issue_df = bootstrap_issue_df,
    bootstrap_issue_top_df = bootstrap_issue_top_df,

    method_key = method_key,

    flag0_df = flag0_df,
    flag1_df = flag1_df,
    flag2_df = flag2_df,
    bootstrap_df = bootstrap_df,
    improper_reason_plot_df = improper_reason_plot_df,
    bootstrap_issue_plot_df = bootstrap_issue_plot_df
  )
}