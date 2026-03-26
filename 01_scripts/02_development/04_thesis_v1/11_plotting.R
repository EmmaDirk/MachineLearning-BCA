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
    "R", "T", "flag",
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
  # label maps
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

  # ------------------------------------------------------------
  # method naming
  # ------------------------------------------------------------
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

  # unique internal id used for grouping, so even if labels collide later,
  # the models are still kept separate
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

  # ------------------------------------------------------------
  # truth lookup
  # ------------------------------------------------------------
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

  df <- results_df %>%
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
    )

  if (drop_flagged) {
    df <- df %>% dplyr::filter(flag == 0)
  }

  df <- df %>%
    dplyr::filter(T %in% occasions)

  # ------------------------------------------------------------
  # long format
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
  # summaries
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
  # ordering
  # ------------------------------------------------------------
  method_key <- df %>%
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
    dplyr::pull(method)

  relbias_df <- relbias_df %>%
    dplyr::mutate(method = factor(method, levels = unique(method_levels)))

  rmse_df <- rmse_df %>%
    dplyr::mutate(method = factor(method, levels = unique(method_levels)))

  se_df <- se_df %>%
    dplyr::mutate(method = factor(method, levels = unique(method_levels)))

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
    # CLPM: saturated lime -> saturated deep blue
    make_family_palette(clpm_methods, "#39FF14", "#0033FF"),

    # RI-CLPM: saturated yellow -> saturated red
    make_family_palette(riclpm_methods, "#FFD600", "#FF1F1F"),

    # DPM: saturated pink -> saturated purple
    make_family_palette(dpm_methods, "#FF1493", "#7A00FF")
  )

  pal_method <- pal_method[unique(method_levels)]

  # ------------------------------------------------------------
  # plotting helper
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
  # plots
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

    relbias_df = relbias_df,
    rmse_df = rmse_df,
    se_df = se_df,
    method_key = method_key,
    pal_method = pal_method,
    plotting_data = long_est_df,
    plotting_se_data = long_se_df
  )
}