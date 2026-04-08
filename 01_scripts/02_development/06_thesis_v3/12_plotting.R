# --------------------------------------------------------------------------------
# Build publication-style plots and summary tables from row-bound simulation output
#
# Returned objects in this version:
#   - summary data frames:
#       relbias_df, rmse_df, se_df, detect_df, se_check,
#       flag0_df, flag1_df, flag2_df, mse_df, r2_df
#   - plots:
#       plot_relbias, plot_rmse, plot_se, plot_power,
#       plot_relbias_clpm, plot_relbias_riclpm, plot_relbias_dpm,
#       plot_rmse_clpm, plot_rmse_riclpm, plot_rmse_dpm,
#       plot_se_clpm, plot_se_riclpm, plot_se_dpm,
#       plot_power_clpm, plot_power_riclpm, plot_power_dpm,
#       plot_se_ratio_clpm, plot_se_ratio_riclpm, plot_se_ratio_dpm,
#       plot_se_diff_clpm, plot_se_diff_riclpm, plot_se_diff_dpm,
#       plot_mse, plot_r2
#
# Main changes relative to the prior version:
#   1. Flag plots are removed from the returned object.
#   2. Flag summary data frames are still returned unchanged.
#   3. SE diagnostic plots are now family-specific instead of combined.
# --------------------------------------------------------------------------------

plot_engine_results <- function(results_df,
                                drop_flagged = TRUE,
                                occasions = 2:5,
                                alpha = 0.05) {

  # ==============================================================================
  # 1. Basic input checks
  # ==============================================================================

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

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).")
  }

  # ==============================================================================
  # 2. Small helper objects and helper functions
  # ==============================================================================

  # Compact model code -> readable model family name
  model_map <- c(
    "C" = "CLPM",
    "R" = "RI-CLPM",
    "D" = "DPM"
  )

  # Compact residualiser code -> readable residualiser name
  resid_map <- c(
    "N" = "None",
    "L" = "LM",
    "X" = "XGB",
    "E" = "EN"
  )

  # Clean helper:
  # turn "" into NA, keep everything as character
  clean_text <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
  }

  # Shared grouped mean + MCSE helper
  # Used for summaries where the target is just mean(value) and MCSE of that mean.
  summarise_mean_mcse <- function(data, value_col, mean_name, mcse_name, group_cols) {
    data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::summarise(
        nsim = dplyr::n(),
        !!mean_name := mean(.data[[value_col]], na.rm = TRUE),
        !!mcse_name := stats::sd(.data[[value_col]], na.rm = TRUE) / sqrt(nsim),
        .groups = "drop"
      )
  }

  # ------------------------------------------------------------------------------
  # Generic line plot helper
  #
  # This is the central plotting engine used by the main/family/diagnostic plots.
  # Arguments:
  #   data        : plotting data frame
  #   y           : name of y column
  #   color_var   : grouping/color variable
  #   group_var   : line grouping variable
  #   facet_var   : optional facet variable
  #   y_se        : optional standard error / MCSE column for error bars
  #   palette     : named vector of colours
  #   y_label     : y axis label
  #   ref_line    : optional horizontal reference line
  #   zero_line   : add y = 0 dashed line
  #   percent_y   : format y axis as percentages
  #   clamp_01    : clamp errorbars and axis to [0, 1]
  # ------------------------------------------------------------------------------
  make_line_plot <- function(data,
                             y,
                             color_var,
                             group_var,
                             facet_var = NULL,
                             y_se = NULL,
                             palette = NULL,
                             y_label = NULL,
                             ref_line = NULL,
                             zero_line = FALSE,
                             percent_y = FALSE,
                             clamp_01 = FALSE) {

    p <- ggplot2::ggplot(
      data = data,
      mapping = ggplot2::aes(
        x = T,
        y = .data[[y]],
        color = .data[[color_var]],
        group = .data[[group_var]]
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2)

    # Add MCSE / uncertainty bars if requested
    if (!is.null(y_se)) {
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

    # Add faceting if requested
    if (!is.null(facet_var)) {
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), nrow = 1, scales = "fixed")
    }

    # Apply manual palette if given
    if (!is.null(palette)) {
      p <- p + ggplot2::scale_color_manual(values = palette, drop = FALSE)
    }

    # X-axis breaks = all observed occasions in the supplied data
    p <- p + ggplot2::scale_x_continuous(breaks = sort(unique(data$T)))

    # Optional y-axis formatting as percent
    if (percent_y) {
      if (clamp_01) {
        p <- p + ggplot2::scale_y_continuous(
          labels = function(x) scales::percent(x, accuracy = 1),
          limits = c(0, 1)
        )
      } else {
        p <- p + ggplot2::scale_y_continuous(
          labels = function(x) scales::percent(x, accuracy = 1)
        )
      }
    }

    # Optional horizontal reference lines
    if (!is.null(ref_line)) {
      p <- p + ggplot2::geom_hline(
        yintercept = ref_line,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    if (zero_line) {
      p <- p + ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.4
      )
    }

    # Shared theme / labelling
    p +
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
  }

  # ------------------------------------------------------------------------------
  # Family-specific wrapper:
  # keep only methods from one SEM family, then call the generic line plot helper
  # ------------------------------------------------------------------------------
  make_family_plot <- function(data,
                               y,
                               y_se,
                               family_name,
                               method_key,
                               pal_method,
                               y_label = NULL,
                               ref_line = NULL,
                               zero_line = FALSE,
                               percent_y = FALSE,
                               clamp_01 = FALSE) {

    family_methods <- method_key %>%
      dplyr::filter(.data$family == .env$family_name) %>%
      dplyr::pull(method) %>%
      unique()

    family_data <- data %>%
      dplyr::filter(method %in% family_methods) %>%
      dplyr::mutate(method = factor(as.character(method), levels = family_methods))

    family_palette <- pal_method[family_methods]

    make_line_plot(
      data = family_data,
      y = y,
      color_var = "method",
      group_var = "method_id",
      facet_var = "estimand",
      y_se = y_se,
      palette = family_palette,
      y_label = y_label,
      ref_line = ref_line,
      zero_line = zero_line,
      percent_y = percent_y,
      clamp_01 = clamp_01
    )
  }

  # ==============================================================================
  # 3. Prepare the master analysis data set
  # ==============================================================================

  # The original code created labels with row-wise pmap().
  # This version does the same work vectorially, which is much faster on large data.

  df_all <- results_df %>%
    dplyr::mutate(
      # Clean common string-like fields once
      exclusion = clean_text(exclusion),
      c_order = clean_text(c_order),
      free_loadings = clean_text(free_loadings),
      bootstrap_B = clean_text(bootstrap_B),
      improper_reason = clean_text(improper_reason),

      # Expand compact model/residualiser codes
      model_full = dplyr::recode(as.character(model), !!!model_map, .default = as.character(model)),
      resid_full = dplyr::recode(as.character(residualizer), !!!resid_map, .default = as.character(residualizer)),

      # Numeric versions used in label construction
      c_order_num = suppressWarnings(as.numeric(c_order)),
      free_loadings_num = suppressWarnings(as.numeric(free_loadings)),

      # Free-loading variants get renamed
      model_full = dplyr::case_when(
        !is.na(free_loadings_num) & free_loadings_num == 1 & model_full == "RI-CLPM" ~ "fRI-CLPM",
        !is.na(free_loadings_num) & free_loadings_num == 1 & model_full == "DPM" ~ "fDPM",
        TRUE ~ model_full
      ),

      # Unadjusted = no residualisation / no BCA wrapper
      is_unadjusted = is.na(resid_full) | resid_full == "None",

      # Base method label
      base_label = dplyr::case_when(
        is.na(resid_full) ~ model_full,
        resid_full == "None" ~ model_full,
        TRUE ~ paste("BCA", resid_full, model_full)
      ),

      # Detail component shown in brackets
      detail_part = dplyr::case_when(
        # Simple unadjusted RI-CLPM / DPM variants: no bracket detail
        is_unadjusted & model_full %in% c("RI-CLPM", "fRI-CLPM", "DPM", "fDPM") ~ NA_character_,

        # Keep CLPM confounder order visible even when unadjusted
        is_unadjusted & model_full == "CLPM" & !is.na(c_order_num) ~ as.character(c_order_num),
        is_unadjusted & model_full == "CLPM" & is.na(c_order_num) ~ "0",

        # Adjustment cases: show order and exclusion if both exist
        !is.na(c_order_num) & !is.na(exclusion) ~ paste0(c_order_num, ":", exclusion),
        !is.na(c_order_num) ~ as.character(c_order_num),
        !is.na(exclusion) ~ exclusion,
        TRUE ~ NA_character_
      ),

      # Final display label
      method = dplyr::if_else(
        is.na(detail_part),
        base_label,
        paste0(base_label, " (", detail_part, ")")
      ),

      # Stable identifier for grouping methods exactly
      method_id = paste(
        paste0("model=", clean_text(model)),
        paste0("resid=", clean_text(residualizer)),
        paste0("excl=", clean_text(exclusion)),
        paste0("corder=", clean_text(c_order)),
        paste0("free=", clean_text(free_loadings)),
        sep = " | "
      ),

      # ML diagnostics collapse over SEM model, so only residualiser matters there
      ml_method = dplyr::recode(
        as.character(residualizer),
        "L" = "LM",
        "E" = "EN",
        "X" = "XGB",
        "N" = NA_character_,
        .default = as.character(residualizer)
      ),
      ml_method_id = paste0("resid=", clean_text(residualizer))
    ) %>%
    dplyr::filter(T %in% occasions)

  # Main analysis subset:
  # optionally keep only admissible runs
  df <- df_all
  if (drop_flagged) {
    df <- df %>% dplyr::filter(analysis_flag == 0)
  }

  # ==============================================================================
  # 4. Build one long estimate/SE table and reuse it everywhere
  # ==============================================================================

  # This is one of the main speed improvements.
  #
  # We build a long table that contains, for each:
  #   - simulation replicate R
  #   - occasion T
  #   - method
  #   - estimand
  # both:
  #   - estimate
  #   - se_estimate
  # plus the matching true parameter value.
  #
  # Once this exists, relative bias / RMSE / SE / power / SE diagnostics all
  # come from this same object.

  sim_long_all <- df_all %>%
    tidyr::pivot_longer(
      cols = c(ARX, ARY, CXY, CYX, se_ARX, se_ARY, se_CXY, se_CYX),
      names_to = c("se_prefix", "estimand"),
      names_pattern = "(se_)?(ARX|ARY|CXY|CYX)",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      value_type = dplyr::if_else(is.na(se_prefix), "estimate", "se_estimate")
    ) %>%
    dplyr::select(-se_prefix) %>%
    tidyr::pivot_wider(
      names_from = value_type,
      values_from = value
    ) %>%
    dplyr::mutate(
      true = dplyr::case_when(
        estimand == "ARX" ~ beta_x,
        estimand == "ARY" ~ beta_y,
        estimand == "CXY" ~ gamma_xy,
        estimand == "CYX" ~ gamma_yx,
        TRUE ~ NA_real_
      )
    )

  sim_long <- sim_long_all
  if (drop_flagged) {
    sim_long <- sim_long %>% dplyr::filter(analysis_flag == 0)
  }

  # ==============================================================================
  # 5. Classic performance summaries
  # ==============================================================================

  # ------------------------------------------------------------------------------
  # Relative bias
  #
  # rel_bias = (mean estimate - true) / true
  # MCSE derived from the SD of the estimate divided by sqrt(n), then rescaled
  # by |true| to put it on the relative-bias scale.
  # ------------------------------------------------------------------------------
  relbias_df <- sim_long %>%
    dplyr::filter(!is.na(estimate)) %>%
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

  # ------------------------------------------------------------------------------
  # RMSE
  #
  # rmse = sqrt(mean((estimate - true)^2))
  #
  # MCSE of RMSE is computed from the delta method:
  #   if m = mean(sq_err), rmse = sqrt(m)
  #   MCSE(rmse) ≈ sqrt(var(sq_err)/n) / (2 * rmse)
  #
  # This version avoids rowwise() and list-columns.
  # ------------------------------------------------------------------------------
  rmse_df <- sim_long %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(
      sq_err = (estimate - true)^2
    ) %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      mean_sq_err = mean(sq_err, na.rm = TRUE),
      var_sq_err = stats::var(sq_err, na.rm = TRUE),
      rmse = sqrt(mean_sq_err),
      mcse_rmse = dplyr::if_else(
        nsim <= 1 | is.na(rmse) | rmse == 0,
        NA_real_,
        sqrt(var_sq_err / nsim) / (2 * rmse)
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(-mean_sq_err, -var_sq_err)

  # ------------------------------------------------------------------------------
  # Mean reported SE
  # ------------------------------------------------------------------------------
  se_df <- summarise_mean_mcse(
    data = sim_long %>% dplyr::filter(!is.na(se_estimate)),
    value_col = "se_estimate",
    mean_name = "mean_se",
    mcse_name = "mcse_mean_se",
    group_cols = c("T", "estimand", "method_id", "method")
  )

  # ------------------------------------------------------------------------------
  # Detection / power summary
  #
  # Two-sided Wald test at level alpha:
  #   detected = |estimate / se_estimate| > z_(1 - alpha/2)
  #
  # If true effect is 0, the plot is a Type I error plot.
  # Otherwise it is a power plot.
  # ------------------------------------------------------------------------------
  crit_z <- stats::qnorm(1 - alpha / 2)

  detect_df <- sim_long %>%
    dplyr::filter(!is.na(estimate), !is.na(se_estimate), se_estimate > 0) %>%
    dplyr::mutate(
      z_value = estimate / se_estimate,
      detected = abs(z_value) > crit_z
    ) %>%
    dplyr::group_by(T, estimand, method_id, method) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      detect_prob = mean(detected, na.rm = TRUE),
      mcse_detect = sqrt(detect_prob * (1 - detect_prob) / nsim),
      error_type = dplyr::if_else(abs(true) < 1e-12, "Type I error", "Power"),
      .groups = "drop"
    )

  # ==============================================================================
  # 6. SE calibration diagnostics
  # ==============================================================================

  # These are defined on admissible main-analysis runs only.
  # That behaviour is preserved.
  sim_long_used <- sim_long_all %>%
    dplyr::filter(
      analysis_flag == 0,
      !is.na(estimate),
      !is.na(se_estimate)
    )

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

  # ==============================================================================
  # 7. Flag summaries
  # ==============================================================================

  # Flag summaries are still computed and returned as data frames.
  # Only the flag plots have been removed from this version.

  flag_df <- df_all %>%
    tidyr::pivot_longer(
      cols = c(flag0, flag1, flag2),
      names_to = "flag_type",
      values_to = "flag_value"
    ) %>%
    dplyr::group_by(T, method_id, method, flag_type) %>%
    dplyr::summarise(
      prop_flag = mean(flag_value, na.rm = TRUE),
      .groups = "drop"
    )

  flag0_df <- flag_df %>% dplyr::filter(flag_type == "flag0") %>% dplyr::select(-flag_type)
  flag1_df <- flag_df %>% dplyr::filter(flag_type == "flag1") %>% dplyr::select(-flag_type)
  flag2_df <- flag_df %>% dplyr::filter(flag_type == "flag2") %>% dplyr::select(-flag_type)

  # ==============================================================================
  # 8. ML diagnostics
  # ==============================================================================

  # ML diagnostics are attached to the residualiser, not the downstream SEM method.
  # So all ML summaries are grouped by ml_method only.
  df_ml <- df %>%
    dplyr::filter(residualizer != "N", !is.na(ml_method))

  ml_long_df <- df_ml %>%
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

  mse_df <- summarise_mean_mcse(
    data = ml_long_df %>% dplyr::filter(metric == "MSE"),
    value_col = "value",
    mean_name = "mean_mse",
    mcse_name = "mcse_mean_mse",
    group_cols = c("T", "target", "ml_method_id", "ml_method")
  )

  r2_df <- summarise_mean_mcse(
    data = ml_long_df %>% dplyr::filter(metric == "R2"),
    value_col = "value",
    mean_name = "mean_r2",
    mcse_name = "mcse_mean_r2",
    group_cols = c("T", "target", "ml_method_id", "ml_method")
  )

  # ==============================================================================
  # 9. Method ordering and palettes
  # ==============================================================================

  # Method ordering is preserved in a readable family-first way:
  # RI-CLPM methods, then DPM, then CLPM, then anything else.
  method_key <- df_all %>%
    dplyr::distinct(method_id, method) %>%
    dplyr::mutate(
      family = dplyr::case_when(
        grepl("^model=R\\b", method_id) ~ "RI-CLPM",
        grepl("^model=D\\b", method_id) ~ "DPM",
        grepl("^model=C\\b", method_id) ~ "CLPM",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::arrange(family, method)

  method_levels <- method_key$method

  # Method palette for SEM methods
  pal_method <- setNames(
    scales::hue_pal()(length(method_levels)),
    method_levels
  )

  # ML residualiser palette
  ml_method_levels <- c("LM", "EN", "XGB")
  pal_ml <- setNames(
    scales::hue_pal()(length(ml_method_levels)),
    ml_method_levels
  )

  # Apply common factor ordering to all returned summary frames that use method
  relbias_df$method <- factor(as.character(relbias_df$method), levels = method_levels)
  rmse_df$method <- factor(as.character(rmse_df$method), levels = method_levels)
  se_df$method <- factor(as.character(se_df$method), levels = method_levels)
  detect_df$method <- factor(as.character(detect_df$method), levels = method_levels)
  se_check$method <- factor(as.character(se_check$method), levels = method_levels)
  flag0_df$method <- factor(as.character(flag0_df$method), levels = method_levels)
  flag1_df$method <- factor(as.character(flag1_df$method), levels = method_levels)
  flag2_df$method <- factor(as.character(flag2_df$method), levels = method_levels)

  mse_df$ml_method <- factor(as.character(mse_df$ml_method), levels = ml_method_levels)
  r2_df$ml_method <- factor(as.character(r2_df$ml_method), levels = ml_method_levels)

  # ==============================================================================
  # 10. Main combined plots
  # ==============================================================================

  plot_relbias <- make_line_plot(
    data = relbias_df,
    y = "rel_bias",
    color_var = "method",
    group_var = "method_id",
    facet_var = "estimand",
    y_se = "mcse_rel_bias",
    palette = pal_method,
    y_label = NULL,
    zero_line = TRUE
  )

  plot_rmse <- make_line_plot(
    data = rmse_df,
    y = "rmse",
    color_var = "method",
    group_var = "method_id",
    facet_var = "estimand",
    y_se = "mcse_rmse",
    palette = pal_method,
    y_label = NULL
  )

  plot_se <- make_line_plot(
    data = se_df,
    y = "mean_se",
    color_var = "method",
    group_var = "method_id",
    facet_var = "estimand",
    y_se = "mcse_mean_se",
    palette = pal_method,
    y_label = NULL
  )

  plot_power <- make_line_plot(
    data = detect_df,
    y = "detect_prob",
    color_var = "method",
    group_var = "method_id",
    facet_var = "estimand",
    y_se = "mcse_detect",
    palette = pal_method,
    y_label = paste0("Detection probability (alpha = ", alpha, ")"),
    ref_line = alpha,
    percent_y = TRUE,
    clamp_01 = TRUE
  )

  # ==============================================================================
  # 11. Family-specific classic plots
  # ==============================================================================

  plot_relbias_clpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method,
    zero_line = TRUE
  )

  plot_relbias_riclpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method,
    zero_line = TRUE
  )

  plot_relbias_dpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method,
    zero_line = TRUE
  )

  plot_rmse_clpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_rmse_riclpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_rmse_dpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_se_clpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_se_riclpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_se_dpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method
  )

  plot_power_clpm <- make_family_plot(
    data = detect_df,
    y = "detect_prob",
    y_se = "mcse_detect",
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = paste0("Detection probability (alpha = ", alpha, ")"),
    ref_line = alpha,
    percent_y = TRUE,
    clamp_01 = TRUE
  )

  plot_power_riclpm <- make_family_plot(
    data = detect_df,
    y = "detect_prob",
    y_se = "mcse_detect",
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = paste0("Detection probability (alpha = ", alpha, ")"),
    ref_line = alpha,
    percent_y = TRUE,
    clamp_01 = TRUE
  )

  plot_power_dpm <- make_family_plot(
    data = detect_df,
    y = "detect_prob",
    y_se = "mcse_detect",
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = paste0("Detection probability (alpha = ", alpha, ")"),
    ref_line = alpha,
    percent_y = TRUE,
    clamp_01 = TRUE
  )

  # ==============================================================================
  # 12. SE diagnostic plots
  # ==============================================================================

  # In this version, SE diagnostic plots are family-specific rather than combined.
  # This matches the family-specific structure already used for bias/RMSE/SE/power.

  plot_se_ratio_clpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE / empirical SD",
    ref_line = 1
  )

  plot_se_ratio_riclpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE / empirical SD",
    ref_line = 1
  )

  plot_se_ratio_dpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE / empirical SD",
    ref_line = 1
  )

  plot_se_diff_clpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE - empirical SD",
    ref_line = 0
  )

  plot_se_diff_riclpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "RI-CLPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE - empirical SD",
    ref_line = 0
  )

  plot_se_diff_dpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "DPM",
    method_key = method_key,
    pal_method = pal_method,
    y_label = "Mean reported SE - empirical SD",
    ref_line = 0
  )

  # ==============================================================================
  # 13. ML plots
  # ==============================================================================

  plot_mse <- make_line_plot(
    data = mse_df,
    y = "mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    facet_var = "target",
    y_se = "mcse_mean_mse",
    palette = pal_ml,
    y_label = "Mean OOF MSE"
  )

  plot_r2 <- make_line_plot(
    data = r2_df,
    y = "mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    facet_var = "target",
    y_se = "mcse_mean_r2",
    palette = pal_ml,
    y_label = "Mean OOF R^2",
    ref_line = 0
  )

  # ==============================================================================
  # 14. Return object
  # ==============================================================================

  # Keep the returned object names unchanged where possible.
  # The only structural change is that:
  #   - flag plots are omitted
  #   - SE diagnostic plots are returned in family-specific form
  list(
    # --------------------------------------------------------------------------
    # Summary data frames
    # --------------------------------------------------------------------------
    relbias_df = relbias_df,
    rmse_df = rmse_df,
    se_df = se_df,
    detect_df = detect_df,
    se_check = se_check,
    flag0_df = flag0_df,
    flag1_df = flag1_df,
    flag2_df = flag2_df,
    mse_df = mse_df,
    r2_df = r2_df,

    # --------------------------------------------------------------------------
    # Classic combined plots
    # --------------------------------------------------------------------------
    plot_relbias = plot_relbias,
    plot_rmse = plot_rmse,
    plot_se = plot_se,
    plot_power = plot_power,

    # --------------------------------------------------------------------------
    # Family-specific classic plots
    # --------------------------------------------------------------------------
    plot_relbias_clpm = plot_relbias_clpm,
    plot_relbias_riclpm = plot_relbias_riclpm,
    plot_relbias_dpm = plot_relbias_dpm,
    plot_rmse_clpm = plot_rmse_clpm,
    plot_rmse_riclpm = plot_rmse_riclpm,
    plot_rmse_dpm = plot_rmse_dpm,
    plot_se_clpm = plot_se_clpm,
    plot_se_riclpm = plot_se_riclpm,
    plot_se_dpm = plot_se_dpm,
    plot_power_clpm = plot_power_clpm,
    plot_power_riclpm = plot_power_riclpm,
    plot_power_dpm = plot_power_dpm,

    # --------------------------------------------------------------------------
    # Family-specific SE diagnostic plots
    # --------------------------------------------------------------------------
    plot_se_ratio_clpm = plot_se_ratio_clpm,
    plot_se_ratio_riclpm = plot_se_ratio_riclpm,
    plot_se_ratio_dpm = plot_se_ratio_dpm,
    plot_se_diff_clpm = plot_se_diff_clpm,
    plot_se_diff_riclpm = plot_se_diff_riclpm,
    plot_se_diff_dpm = plot_se_diff_dpm,

    # --------------------------------------------------------------------------
    # ML plots
    # --------------------------------------------------------------------------
    plot_mse = plot_mse,
    plot_r2 = plot_r2
  )
}