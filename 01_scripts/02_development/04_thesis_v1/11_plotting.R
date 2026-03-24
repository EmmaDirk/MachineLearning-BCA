# This script takes the row-bound simulation output and builds publication-style plots
# for:
# - relative bias
# - RMSE
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
    "ARX", "ARY", "CXY", "CYX"
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
    bootstrap_num <- suppressWarnings(as.numeric(bootstrap_B))

    base_label <- dplyr::case_when(
      is.na(resid_full) ~ model_full,
      resid_full == "None" ~ model_full,
      TRUE ~ paste("BCA", resid_full, model_full)
    )

    # IMPORTANT:
    # always print the confounder order so that order 0 and order 1
    # never collapse onto the same plotted label
    order_part <- dplyr::case_when(
      is.na(c_order_num) ~ NA_character_,
      TRUE ~ paste0("order ", c_order_num)
    )

    loading_part <- dplyr::case_when(
      is.na(free_loadings_num) ~ NA_character_,
      free_loadings_num == 1 ~ "free loadings",
      free_loadings_num == 0 ~ NA_character_,
      TRUE ~ NA_character_
    )

    bootstrap_part <- dplyr::case_when(
      is.na(bootstrap_num) ~ NA_character_,
      bootstrap_num >= 2 ~ paste0("boot ", bootstrap_num),
      TRUE ~ NA_character_
    )

    exclusion_part <- dplyr::case_when(
      is.na(exclusion) ~ NA_character_,
      TRUE ~ paste0("excl: ", exclusion)
    )

    extras <- c(order_part, loading_part, bootstrap_part, exclusion_part)
    extras <- extras[!is.na(extras)]

    if (length(extras) == 0) {
      base_label
    } else {
      paste0(base_label, " (", paste(extras, collapse = ", "), ")")
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
      paste0("boot=", clean_text(bootstrap_B)),
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
  long_df <- df %>%
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

  # ------------------------------------------------------------
  # summaries
  # ------------------------------------------------------------
  relbias_df <- long_df %>%
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

  rmse_df <- long_df %>%
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

  # ------------------------------------------------------------
  # ordering
  # ------------------------------------------------------------
  method_levels <- relbias_df %>%
    dplyr::distinct(method_id, method) %>%
    dplyr::arrange(method) %>%
    dplyr::pull(method)

  relbias_df <- relbias_df %>%
    dplyr::mutate(method = factor(method, levels = unique(method_levels)))

  rmse_df <- rmse_df %>%
    dplyr::mutate(method = factor(method, levels = unique(method_levels)))

  # ------------------------------------------------------------
  # palette
  # ------------------------------------------------------------
  pal_method <- viridisLite::viridis(
    n = length(unique(method_levels)),
    option = "viridis",
    begin = 0.1,
    end = 0.9
  )
  names(pal_method) <- unique(method_levels)

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
      ggplot2::scale_color_manual(values = pal_method) +
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

  # ------------------------------------------------------------
  # method key
  # ------------------------------------------------------------
  method_key <- df %>%
    dplyr::distinct(
      model, residualizer, exclusion, c_order,
      free_loadings, bootstrap_B, method_id, method
    ) %>%
    dplyr::arrange(method)

  # ------------------------------------------------------------
  # return
  # ------------------------------------------------------------
  list(
    plot_relbias = plot_relbias,
    plot_rmse = plot_rmse,
    relbias_df = relbias_df,
    rmse_df = rmse_df,
    method_key = method_key,
    pal_method = pal_method,
    plotting_data = long_df
  )
}