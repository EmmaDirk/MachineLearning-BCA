# ============================================================
# Adaptive plotting for the new simulation engine output
# ============================================================

plot_engine_results <- function(results_df,
                                estimand_to_plot = "CXY",
                                drop_flagged = TRUE,
                                scenario_filter = NULL) {

  # ------------------------------------------------------------
  # checks
  # ------------------------------------------------------------
  required_cols <- c(
    "R", "T", "flag",
    "model", "residualizer", "exclusion", "c_order",
    "beta_x", "beta_y", "gamma_xy", "gamma_yx",
    "ARX", "ARY", "CXY", "CYX"
  )

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ------------------------------------------------------------
  # helper maps
  # ------------------------------------------------------------
  model_map <- c(
    "C" = "CLPM",
    "R" = "RI-CLPM",
    "D" = "DPM"
  )

  resid_map <- c(
    "N"    = "None",
    "none" = "None",
    "L"    = "LM",
    "linear" = "LM",
    "X"    = "XGB",
    "xgb"  = "XGB"
  )

  clean_text <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
  }

  # ------------------------------------------------------------
  # helper: construct readable method label
  # ------------------------------------------------------------
  make_method_label <- function(model, residualizer, exclusion, c_order) {

    model_full <- dplyr::recode(as.character(model), !!!model_map, .default = as.character(model))
    resid_full <- dplyr::recode(as.character(residualizer), !!!resid_map, .default = as.character(residualizer))

    exclusion <- clean_text(exclusion)
    c_order   <- clean_text(c_order)

    # baseline covariate adjustment part
    bca_part <- dplyr::case_when(
      is.na(resid_full) ~ NA_character_,
      resid_full %in% c("None", "N") ~ NA_character_,
      TRUE ~ paste("BCA", resid_full)
    )

    # base label
    base_label <- dplyr::case_when(
      is.na(bca_part) ~ model_full,
      TRUE            ~ paste(bca_part, model_full)
    )

    # order label:
    # keep it generic, because in your setup the interpretation of c_order
    # depends on the model / residualizer combination
    order_part <- dplyr::case_when(
      is.na(c_order) ~ NA_character_,
      TRUE ~ paste0("order ", c_order)
    )

    exclusion_part <- dplyr::case_when(
      is.na(exclusion) ~ NA_character_,
      TRUE ~ paste0("excl: ", exclusion)
    )

    extras <- c(order_part, exclusion_part)
    extras <- extras[!is.na(extras)]

    if (length(extras) == 0) {
      base_label
    } else {
      paste0(base_label, " (", paste(extras, collapse = ", "), ")")
    }
  }

  # ------------------------------------------------------------
  # helper: true parameter lookup
  # ------------------------------------------------------------
  truth_map <- tibble::tibble(
    estimand = c("ARX", "ARY", "CXY", "CYX"),
    true_col = c("beta_x", "beta_y", "gamma_xy", "gamma_yx")
  )

  # ------------------------------------------------------------
  # prepare data
  # ------------------------------------------------------------
  df <- results_df %>%
    dplyr::mutate(
      exclusion   = clean_text(exclusion),
      c_order     = clean_text(c_order),
      model_label = dplyr::recode(as.character(model), !!!model_map, .default = as.character(model)),
      resid_label = dplyr::recode(as.character(residualizer), !!!resid_map, .default = as.character(residualizer))
    )

  if (drop_flagged) {
    df <- df %>% dplyr::filter(.data$flag == 0)
  }

  # Each distinct true setup becomes a scenario.
  # T is included so runs with different total wave counts do not get mixed.
  df <- df %>%
    dplyr::mutate(
      scenario_id = dplyr::dense_rank(
        dplyr::pick(T, beta_x, beta_y, gamma_xy, gamma_yx)
      ),
      scenario_label = paste0(
        "Scenario ", scenario_id,
        "\nT=", T,
        ", βx=", beta_x,
        ", βy=", beta_y,
        ", γxy=", gamma_xy,
        ", γyx=", gamma_yx
      ),
      method = purrr::pmap_chr(
        list(model, residualizer, exclusion, c_order),
        make_method_label
      )
    )

  # optional scenario filter
  if (!is.null(scenario_filter)) {
    if (is.numeric(scenario_filter)) {
      df <- df %>% dplyr::filter(.data$scenario_id %in% scenario_filter)
    } else {
      df <- df %>% dplyr::filter(.data$scenario_label %in% scenario_filter)
    }
  }

  # ------------------------------------------------------------
  # long format for estimates
  # ------------------------------------------------------------
  long_df <- df %>%
    dplyr::filter(.data$T >= 2) %>%
    tidyr::pivot_longer(
      cols = c("ARX", "ARY", "CXY", "CYX"),
      names_to = "estimand",
      values_to = "estimate"
    ) %>%
    dplyr::filter(!is.na(.data$estimate)) %>%
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

  # ------------------------------------------------------------
  # summaries
  # ------------------------------------------------------------
  relbias_df <- long_df %>%
    dplyr::group_by(scenario_id, scenario_label, T, estimand, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      mean_est = mean(estimate, na.rm = TRUE),
      rel_bias = (mean_est - true) / true,
      mcse_rel_bias = stats::sd(estimate, na.rm = TRUE) / sqrt(nsim) / abs(true),
      .groups = "drop"
    )

  rmse_df <- long_df %>%
    dplyr::group_by(scenario_id, scenario_label, T, estimand, method) %>%
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

  # ------------------------------------------------------------
  # consistent method ordering
  # ------------------------------------------------------------
  method_levels <- relbias_df %>%
    dplyr::distinct(method) %>%
    dplyr::arrange(method) %>%
    dplyr::pull(method)

  relbias_df <- relbias_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  rmse_df <- rmse_df %>%
    dplyr::mutate(method = factor(method, levels = method_levels))

  # ------------------------------------------------------------
  # palette
  # ------------------------------------------------------------
  pal_method <- viridisLite::viridis(
    n = length(method_levels),
    option = "viridis",
    begin = 0.1,
    end = 0.9
  )
  names(pal_method) <- method_levels

  # ------------------------------------------------------------
  # plotting helper
  # ------------------------------------------------------------
  make_metric_plot <- function(data, metric, metric_se, ylab, add_zero_line = FALSE) {

    p <- data %>%
      dplyr::filter(.data$estimand == estimand_to_plot) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = T,
          y = .data[[metric]],
          color = method,
          group = method
        )
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data[[metric]] - .data[[metric_se]],
          ymax = .data[[metric]] + .data[[metric_se]]
        ),
        width = 0.15,
        linewidth = 0.4
      ) +
      ggplot2::scale_color_manual(values = pal_method) +
      ggplot2::labs(
        x = "Occasion",
        y = ylab,
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        panel.spacing.x = grid::unit(1.2, "cm"),
        strip.text = ggplot2::element_text(size = 13),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )

    if (add_zero_line) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4)
    }

    # facet by the true generating setup
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_wrap2(~ scenario_label, nrow = 1, axes = "y")
    } else {
      p <- p + ggplot2::facet_wrap(~ scenario_label, nrow = 1, scales = "free_y")
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
    ylab = "Relative bias",
    add_zero_line = TRUE
  )

  plot_rmse <- make_metric_plot(
    data = rmse_df,
    metric = "rmse",
    metric_se = "mcse_rmse",
    ylab = "RMSE",
    add_zero_line = FALSE
  )

  plot_relbias_noleg <- plot_relbias +
    ggplot2::theme(legend.position = "none")

  combined_plot <- plot_relbias_noleg / plot_rmse

  # ------------------------------------------------------------
  # method key table
  # ------------------------------------------------------------
  method_key <- df %>%
    dplyr::distinct(
      model, residualizer, exclusion, c_order, method
    ) %>%
    dplyr::arrange(method)

  # ------------------------------------------------------------
  # return
  # ------------------------------------------------------------
  list(
    combined_plot = combined_plot,
    plot_relbias = plot_relbias,
    plot_rmse = plot_rmse,
    relbias_df = relbias_df,
    rmse_df = rmse_df,
    method_key = method_key,
    pal_method = pal_method,
    plotting_data = long_df
  )
}