# =================================================================================================
#
# These helpers create reusable overviews for simulation-results data:
# - overview_models() gives a compact table of distinct model specifications
# - plot_overview_suite() creates diagnostic summary tables and plots
#
# The executable R code has not been changed. Only comments have been reorganized.
# =================================================================================================

# ---- overview_models -----------------------------------------------------------------------------

overview_models <- function(
    results_df,
    model_cols = c(
      "model",
      "residualizer",
      "sem_exclusion",
      "sem_c_order",
      "residualizer_exclusion",
      "residualizer_c_order",
      "free_loadings"
    ),
    include_counts = FALSE,
    sort_by = c(
      "model",
      "residualizer",
      "sem_exclusion",
      "sem_c_order",
      "residualizer_exclusion",
      "residualizer_c_order",
      "free_loadings"
    ),
    print_result = TRUE
) {
# ---- purpose -------------------------------------------------------------------------------------
  # This function gives a compact overview of which distinct model setups are
  # present in a simulation-results data frame.
  #
  # It returns one row per unique model specification, based only on the model
  # columns you care about.
  #
  # This version is intentionally stripped down:
  # - no bootstrap_B
  # - no model_id
  # - no model_name
  # - no model_logic
  # - no count columns such as n_rows, n_replications, T_min, T_max
  #
  # So the output is just a clean table of unique model-defining combinations.

# ---- 1. basic input checks -----------------------------------------------------------------------
  stopifnot(is.data.frame(results_df))

  missing_model_cols <- setdiff(model_cols, names(results_df))
  if (length(missing_model_cols) > 0) {
    stop(
      "These model-defining columns are missing from results_df: ",
      paste(missing_model_cols, collapse = ", ")
    )
  }

# ---- 2. small helper functions -------------------------------------------------------------------
  clean_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
  }

  # Convert internal SEM model codes into readable labels.
  model_label <- function(model_code, free_loadings = NA) {
    model_code <- clean_chr(model_code)

    dplyr::case_when(
      model_code == "C" & !is.na(free_loadings) & free_loadings == 1 ~ "fCLPM",
      model_code == "C" ~ "CLPM",
      model_code == "R" & !is.na(free_loadings) & free_loadings == 1 ~ "fRI-CLPM",
      model_code == "R" ~ "RI-CLPM",
      model_code == "D" & !is.na(free_loadings) & free_loadings == 1 ~ "fDPM",
      model_code == "D" ~ "DPM",
      TRUE ~ model_code
    )
  }

  # Convert residualizer codes into readable labels.
  resid_label <- function(resid_code) {
    resid_code <- clean_chr(resid_code)

    dplyr::case_when(
      resid_code == "N" ~ "None",
      resid_code == "L" ~ "LM",
      resid_code == "E" ~ "EN",
      resid_code == "X" ~ "XGB",
      TRUE ~ resid_code
    )
  }

  # Format exclusion fields for display.
  # Missing values are shown as "no".
  fmt_excl <- function(x) {
    x <- clean_chr(x)
    ifelse(is.na(x), "no", x)
  }

  # Convert order variables to character for clean printing.
  fmt_order <- function(x) {
    ifelse(is.na(x), NA_character_, as.character(x))
  }

# ---- 3. build the table of unique model combinations ---------------------------------------------
  # Reduce the full results data frame to one row per distinct model setup.
  model_df <- results_df %>%
    dplyr::mutate(
      model = clean_chr(.data$model),
      residualizer = clean_chr(.data$residualizer),
      sem_exclusion = clean_chr(.data$sem_exclusion),
      residualizer_exclusion = clean_chr(.data$residualizer_exclusion)
    ) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(model_cols)))

# ---- 4. add human-readable labels ----------------------------------------------------------------
  # The raw columns are useful, but the overview is easier to inspect if we also
  # add readable labels for the model family and residualizer.
  model_df <- model_df %>%
    dplyr::mutate(
      model_tag = model_label(model, free_loadings),
      residualizer_tag = resid_label(residualizer),
      sem_exclude = fmt_excl(sem_exclusion),
      resid_exclude = fmt_excl(residualizer_exclusion),
      sem_order = fmt_order(sem_c_order),
      resid_order = fmt_order(residualizer_c_order),
      free_loadings = ifelse(is.na(free_loadings), NA, free_loadings)
    )

# ---- 5. optional count block ---------------------------------------------------------------------
  # Kept here only so the function structure stays flexible, but by default this
  # version does not add count summaries.
  if (include_counts) {
    warning("include_counts = TRUE is currently ignored in this stripped-down version.")
  }

# ---- 6. order and clean the final output table ---------------------------------------------------
  preferred_order <- c(
    "model_tag",
    "residualizer_tag",
    "sem_exclude",
    "sem_order",
    "resid_exclude",
    "resid_order",
    "free_loadings"
  )

  keep_cols <- intersect(preferred_order, names(model_df))

  sort_by <- intersect(sort_by, names(model_df))
  if (length(sort_by) > 0) {
    model_df <- model_df %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(sort_by)))
  }

  model_df <- model_df %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    tibble::as_tibble()

# ---- 7. print and return -------------------------------------------------------------------------
  if (print_result) {
    print(model_df, n = Inf, width = Inf)
  }

  invisible(model_df)
}


# ---- plot_overview_suite -------------------------------------------------------------------------

plot_overview_suite <- function(
    results_df,
    drop_flagged = TRUE,
    occasions = NULL,
    alpha = 0.05,
    improper_top_n = 5,
    print_messages = TRUE
) {
# ---- purpose -------------------------------------------------------------------------------------
  # This function builds a fast overview-plot suite for one simulation-results
  # data frame.
  #
  # It is intended for repeated use across many result files, so the structure is:
  #
  #   1. validate and lightly standardise the input once,
  #   2. build a few reusable long-format tables once,
  #   3. compute all summary tables from those cached tables,
  #   4. create the overview plots from the summaries.
  #
  # The goal is speed and readability. Expensive wrangling is done once and then
  # reused, rather than repeating reshaping work for each plot.
  #
# ---- main outputs --------------------------------------------------------------------------------
  # Returned data frames:
  #   - relbias_df
  #   - rmse_df
  #   - se_df
  #   - detect_df
  #   - se_check
  #   - mse_df
  #   - r2_df
  #   - theory_df
  #   - true_delta_df
  #   - flag0_df
  #   - flag1_df
  #   - flag2_df
  #   - improper_reason_df
  #
  # Returned plots:
  #   - family-specific relative-bias plots
  #   - family-specific RMSE plots
  #   - family-specific mean-SE plots
  #   - family-specific power plots
  #   - family-specific SE-calibration plots (ratio and difference)
  #   - MSE plots for X and Y with theoretical lower bound
  #   - R² plots for X and Y with true-R² reference line
  #   - true-delta plot
  #   - horizontal bar plots for flag0 / flag1 / flag2
  #   - faceted improper-reason plot
  #
# ---- arguments -----------------------------------------------------------------------------------
  # results_df:
  #   Simulation-results data frame.
  #
  # drop_flagged:
  #   If TRUE, rows with analysis_flag != 0 are excluded from the main parameter
  #   and ML summaries. Flag summaries themselves always use the full run-level
  #   data because that is the whole point of those plots.
  #
  # occasions:
  #   Optional numeric vector of T values to keep. If NULL, all observed T values
  #   are kept. In many cases you may want 2:5 because T = 1 often has NA values
  #   for lagged-effect estimates.
  #
  # alpha:
  #   Significance level used for the detection / power summaries.
  #
  # improper_top_n:
  #   Number of improper-solution reasons to keep per method in the improper-
  #   reason plot.
  #
  # print_messages:
  #   If TRUE, print a few light diagnostic messages and warnings.
  #
# ---- 1. basic input checks -----------------------------------------------------------------------
  stopifnot(is.data.frame(results_df))

  required_cols <- c(
    "R", "T", "analysis_flag",
    "flag0", "flag1", "flag2",
    "improper_reason",
    "model", "residualizer",
    "sem_exclusion", "sem_c_order",
    "residualizer_exclusion", "residualizer_c_order",
    "free_loadings",
    "beta_x", "beta_y", "gamma_xy", "gamma_yx",
    "true_r2_x", "true_r2_y",
    "mse_x", "mse_y", "r2_x", "r2_y",
    "ARX", "ARY", "CXY", "CYX",
    "se_ARX", "se_ARY", "se_CXY", "se_CYX"
  )

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.null(occasions) && !is.numeric(occasions)) {
    stop("occasions must be NULL or a numeric vector.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).")
  }

  if (!is.numeric(improper_top_n) || length(improper_top_n) != 1 || improper_top_n < 1) {
    stop("improper_top_n must be a single positive number.")
  }

  # If multiple scenarios or sample sizes are present, this overview will pool
  # across them. That may be intended, but it is worth warning about.
  scenario_cols_present <- intersect(c("scenario_id", "scenario_label", "N"), names(results_df))
  if (length(scenario_cols_present) > 0 && print_messages) {
    varying_counts <- vapply(
      scenario_cols_present,
      function(v) dplyr::n_distinct(results_df[[v]], na.rm = TRUE),
      numeric(1)
    )

    if (any(varying_counts > 1)) {
      message(
        "Note: this data frame contains multiple values for: ",
        paste(names(varying_counts)[varying_counts > 1], collapse = ", "),
        ". This overview aggregates across them."
      )
    }
  }

# ---- 2. small helper functions -------------------------------------------------------------------
  # These helpers keep the main body readable and ensure consistent logic across
  # all summaries and plots.

  # Convert input to plain character and treat empty strings as missing.
  clean_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
  }

  # Convert compact SEM model codes to readable labels.
  model_label <- function(model_code, free_loadings = NA) {
    model_code <- clean_chr(model_code)

    dplyr::case_when(
      model_code == "C" & !is.na(free_loadings) & free_loadings == 1 ~ "fCLPM",
      model_code == "C" ~ "CLPM",
      model_code == "R" & !is.na(free_loadings) & free_loadings == 1 ~ "fRI-CLPM",
      model_code == "R" ~ "RI-CLPM",
      model_code == "D" & !is.na(free_loadings) & free_loadings == 1 ~ "fDPM",
      model_code == "D" ~ "DPM",
      TRUE ~ model_code
    )
  }

  # Convert residualizer codes to readable labels.
  resid_label <- function(resid_code) {
    resid_code <- clean_chr(resid_code)

    dplyr::case_when(
      resid_code == "N" ~ "None",
      resid_code == "L" ~ "LM",
      resid_code == "E" ~ "EN",
      resid_code == "X" ~ "XGB",
      TRUE ~ resid_code
    )
  }

  # Missing exclusion values are displayed as "no".
  fmt_excl <- function(x) {
    x <- clean_chr(x)
    ifelse(is.na(x), "no", x)
  }

  # Convert order values to character for display use.
  fmt_order <- function(x) {
    ifelse(is.na(x), NA_character_, as.character(x))
  }

  # Generic grouped mean + MCSE helper:
  #   mean(V)
  #   MCSE(mean(V)) = sd(V) / sqrt(n)
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

  # Helpers for ordering labels separately within facets in the improper-reason
  # plot, without introducing another package dependency.
  reorder_within_base <- function(x, by, within, fun = mean, sep = "___") {
    stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
  }

  strip_reorder_suffix <- function(x, sep = "___") {
    gsub(paste0(sep, ".*$"), "", x)
  }

# --------------------------------------------------------------------------------------------------
  # Generic line-plot helper used across the overview suite.
  #
  # Features:
  # - optional uncertainty bars,
  # - optional constant reference line,
  # - optional time-varying reference line from a separate data frame,
  # - shared visual style across the plot suite.
# --------------------------------------------------------------------------------------------------
  make_line_plot <- function(data,
                             y,
                             color_var,
                             group_var,
                             facet_var = NULL,
                             y_se = NULL,
                             palette = NULL,
                             y_label = NULL,
                             ref_line = NULL,
                             ref_df = NULL,
                             ref_y = NULL,
                             ref_linetype = "dotted",
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

    # Optional uncertainty bars.
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

    # Optional time-varying theoretical / true reference line.
    if (!is.null(ref_df) && !is.null(ref_y)) {
      p <- p +
        ggplot2::geom_line(
          data = ref_df,
          mapping = ggplot2::aes(x = T, y = .data[[ref_y]]),
          inherit.aes = FALSE,
          color = "black",
          linetype = ref_linetype,
          linewidth = 0.7
        )
    }

    # Optional constant reference line.
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

    if (!is.null(facet_var)) {
      p <- p +
        ggplot2::facet_wrap(
          stats::as.formula(paste("~", facet_var)),
          nrow = 1,
          scales = "fixed"
        )
    }

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
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )
  }

  # Family-specific wrapper around the generic line-plot helper.
  make_family_plot <- function(data,
                               y,
                               family_name,
                               method_key,
                               palette,
                               y_se = NULL,
                               y_label = NULL,
                               ref_line = NULL,
                               zero_line = FALSE,
                               percent_y = FALSE,
                               clamp_01 = FALSE) {

    family_methods <- method_key %>%
      dplyr::filter(family == .env$family_name) %>%
      dplyr::pull(method) %>%
      unique()

    family_data <- data %>%
      dplyr::filter(method %in% family_methods) %>%
      dplyr::mutate(method = factor(as.character(method), levels = family_methods))

    family_palette <- palette[family_methods]

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

  # Generic horizontal bar plot helper for the flag summaries.
  make_flag_bar_plot <- function(data, x, x_label, decreasing = TRUE, palette = NULL) {
    plot_df <- data %>%
      dplyr::arrange(if (decreasing) dplyr::desc(.data[[x]]) else .data[[x]]) %>%
      dplyr::mutate(method = factor(method, levels = method))

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[x]],
        y = method,
        fill = method
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
        axis.text = ggplot2::element_text(size = 11)
      )

    if (!is.null(palette)) {
      p <- p + ggplot2::scale_fill_manual(values = palette, drop = FALSE)
    }

    p
  }

# --------------------------------------------------------------------------------------------------
  # Helper to extract the true delta trajectories from the list column
  # true_delta_t_vector.
  #
  # This is a purely "truth" based object. It does not depend on fitted model
  # performance at all. It only reshapes the stored data-generating deltas into
  # a tidy table that can be plotted.
  #
  # Expected input structure:
  #   each row contains a named numeric vector such as
  #   c("x__c1" = ..., "x__c2" = ..., "y__c1" = ..., ...)
  #
  # Output columns:
  #   - T
  #   - delta_name
  #   - delta
  #   - outcome      (e.g. "x" or "y")
  #   - confounder   (e.g. "c1", "c2", ...)
  #
  # Distinct rows are used first because the same true delta vector is often
  # repeated many times in the simulation output.
# --------------------------------------------------------------------------------------------------
  extract_true_delta_df <- function(data) {
    if (!"true_delta_t_vector" %in% names(data)) {
      return(tibble::tibble())
    }

    delta_base <- data %>%
      dplyr::select(T, true_delta_t_vector) %>%
      dplyr::filter(!purrr::map_lgl(true_delta_t_vector, is.null)) %>%
      dplyr::distinct()

    if (nrow(delta_base) == 0) {
      return(tibble::tibble())
    }

    delta_df <- delta_base %>%
      tidyr::unnest_longer(
        true_delta_t_vector,
        values_to = "delta",
        indices_to = "delta_name"
      ) %>%
      tidyr::separate(
        delta_name,
        into = c("outcome", "confounder"),
        sep = "__",
        remove = FALSE,
        fill = "right"
      ) %>%
      dplyr::mutate(
        outcome = clean_chr(outcome),
        confounder = clean_chr(confounder)
      ) %>%
      dplyr::arrange(outcome, confounder, T)

    delta_df
  }

# ---- 3. standardise / prepare the master table once ----------------------------------------------
  # This is the central wrangling block.
  #
  # The idea is to create the reusable display columns and IDs once and then
  # reuse them across all summaries and plots.

  df_all <- results_df %>%
    dplyr::mutate(
      model = clean_chr(model),
      residualizer = clean_chr(residualizer),
      sem_exclusion = clean_chr(sem_exclusion),
      residualizer_exclusion = clean_chr(residualizer_exclusion),
      improper_reason = clean_chr(improper_reason),

      model_tag = model_label(model, free_loadings),
      residualizer_tag = resid_label(residualizer),

      family = dplyr::case_when(
        model == "C" ~ "CLPM",
        model == "R" ~ "RI-CLPM",
        model == "D" ~ "DPM",
        TRUE ~ "Other"
      ),

      sem_excl_disp = fmt_excl(sem_exclusion),
      resid_excl_disp = fmt_excl(residualizer_exclusion),
      sem_order_disp = fmt_order(sem_c_order),
      resid_order_disp = fmt_order(residualizer_c_order),

      # Compact readable method label used in legends and bar plots.
      method = dplyr::case_when(
        residualizer_tag == "None" ~
          paste0(
            model_tag,
            " [sem ", dplyr::coalesce(sem_order_disp, "NA"),
            "; excl=", sem_excl_disp, "]"
          ),
        TRUE ~
          paste0(
            residualizer_tag,
            " [resid ", dplyr::coalesce(resid_order_disp, "NA"),
            "; excl=", resid_excl_disp, "] + ",
            model_tag,
            " [sem ", dplyr::coalesce(sem_order_disp, "NA"),
            "; excl=", sem_excl_disp, "]"
          )
      ),

      # Stable exact-method identifier.
      method_id = paste(
        paste0("model=", clean_chr(model)),
        paste0("resid=", clean_chr(residualizer)),
        paste0("sem_excl=", clean_chr(sem_exclusion)),
        paste0("sem_c=", sem_c_order),
        paste0("resid_excl=", clean_chr(residualizer_exclusion)),
        paste0("resid_c=", residualizer_c_order),
        paste0("free=", free_loadings),
        sep = " | "
      ),

      # ML-only method label and ID.
      ml_method = residualizer_tag,
      ml_method_id = paste0("resid=", clean_chr(residualizer))
    )

  if (!is.null(occasions)) {
    df_all <- df_all %>% dplyr::filter(T %in% occasions)
  }

  # Main subset used for parameter and ML summaries.
  # Flag summaries use the separate run-level table below.
  df_main <- df_all
  if (drop_flagged) {
    df_main <- df_main %>% dplyr::filter(analysis_flag == 0)
  }

# ---- 4. build the long parameter table once ------------------------------------------------------
  # Backbone for:
  # - relative bias,
  # - RMSE,
  # - mean SE,
  # - power / detection,
  # - SE calibration.

  df_all_id <- df_all %>%
    dplyr::mutate(.row_id = dplyr::row_number())

  est_long <- df_all_id %>%
    tidyr::pivot_longer(
      cols = c(ARX, ARY, CXY, CYX),
      names_to = "estimand",
      values_to = "estimate"
    )

  se_long <- df_all_id %>%
    tidyr::pivot_longer(
      cols = c(se_ARX, se_ARY, se_CXY, se_CYX),
      names_to = "estimand",
      names_prefix = "se_",
      values_to = "se_estimate"
    ) %>%
    dplyr::select(.row_id, estimand, se_estimate)

  sim_long_all <- est_long %>%
    dplyr::left_join(se_long, by = c(".row_id", "estimand")) %>%
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

# ---- 5. build the ml long table once -------------------------------------------------------------
  # Reusable table for MSE / R² summaries and plots.
  # Only actual residualizer methods are relevant here.

  ml_long <- df_main %>%
    dplyr::filter(residualizer != "N", !is.na(ml_method), ml_method != "None") %>%
    tidyr::pivot_longer(
      cols = c(mse_x, mse_y, r2_x, r2_y),
      names_to = "metric_name",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::case_when(
        metric_name %in% c("mse_x", "mse_y") ~ "MSE",
        metric_name %in% c("r2_x", "r2_y") ~ "R2",
        TRUE ~ NA_character_
      ),
      target = dplyr::case_when(
        metric_name %in% c("mse_x", "r2_x") ~ "X",
        metric_name %in% c("mse_y", "r2_y") ~ "Y",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(metric), !is.na(target), !is.na(value))

# ---- 6. build the run-level table once for flags / improper reasons ------------------------------
  # flag0 / flag1 / flag2 are run-level outcomes, not really time-varying
  # outcomes. The raw data usually repeats them across T, so we collapse to one
  # row per replicate-method combination.

  run_level <- df_all %>%
    dplyr::group_by(method_id, method, family, R) %>%
    dplyr::summarise(
      flag0 = dplyr::first(flag0),
      flag1 = dplyr::first(flag1),
      flag2 = dplyr::first(flag2),
      improper_reason = dplyr::first(improper_reason),
      .groups = "drop"
    )

# ---- 6b. build the true-delta table once ---------------------------------------------------------
  # This table is independent of the fitted model summaries. It only reflects the
  # stored data-generating delta values over time.

  true_delta_df <- extract_true_delta_df(df_all)

  if (nrow(true_delta_df) > 0 && print_messages) {
    n_distinct_delta_vectors <- df_all %>%
      dplyr::select(T, true_delta_t_vector) %>%
      dplyr::distinct() %>%
      nrow()

    if (n_distinct_delta_vectors == 1) {
      message(
        "Note: true_delta_t_vector appears constant across the retained rows. ",
        "The delta plot will therefore show flat lines over time."
      )
    }
  }

# ---- 7. summary tables ---------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
  # 7A. Relative bias
# --------------------------------------------------------------------------------------------------
  relbias_df <- sim_long %>%
    dplyr::filter(!is.na(estimate), !is.na(true)) %>%
    dplyr::group_by(T, estimand, method_id, method, family) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      mean_est = mean(estimate, na.rm = TRUE),
      rel_bias = dplyr::if_else(
        abs(true) < 1e-12,
        NA_real_,
        (mean_est - true) / true
      ),
      mcse_rel_bias = dplyr::if_else(
        abs(true) < 1e-12,
        NA_real_,
        stats::sd(estimate, na.rm = TRUE) / sqrt(nsim) / abs(true)
      ),
      .groups = "drop"
    )

# --------------------------------------------------------------------------------------------------
  # 7B. RMSE
# --------------------------------------------------------------------------------------------------
  rmse_df <- sim_long %>%
    dplyr::filter(!is.na(estimate), !is.na(true)) %>%
    dplyr::mutate(sq_err = (estimate - true)^2) %>%
    dplyr::group_by(T, estimand, method_id, method, family) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
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

# --------------------------------------------------------------------------------------------------
  # 7C. Mean reported SE
# --------------------------------------------------------------------------------------------------
  se_df <- summarise_mean_mcse(
    data = sim_long %>%
      dplyr::filter(!is.na(se_estimate)),
    value_col = "se_estimate",
    mean_name = "mean_se",
    mcse_name = "mcse_mean_se",
    group_cols = c("T", "estimand", "method_id", "method", "family")
  )

# --------------------------------------------------------------------------------------------------
  # 7D. Detection / power summary
  #
  # detected = |estimate / se_estimate| > z_(1 - alpha/2)
  #
  # If true = 0, the result is a Type I error rate.
  # If true != 0, the result is power.
# --------------------------------------------------------------------------------------------------
  crit_z <- stats::qnorm(1 - alpha / 2)

  detect_df <- sim_long %>%
    dplyr::filter(!is.na(estimate), !is.na(se_estimate), se_estimate > 0) %>%
    dplyr::mutate(
      z_value = estimate / se_estimate,
      detected = abs(z_value) > crit_z
    ) %>%
    dplyr::group_by(T, estimand, method_id, method, family) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      true = dplyr::first(true),
      detect_prob = mean(detected, na.rm = TRUE),
      mcse_detect = sqrt(detect_prob * (1 - detect_prob) / nsim),
      error_type = dplyr::if_else(abs(true) < 1e-12, "Type I error", "Power"),
      .groups = "drop"
    )

# --------------------------------------------------------------------------------------------------
  # 7E. SE calibration check
# --------------------------------------------------------------------------------------------------
  se_check <- sim_long %>%
    dplyr::filter(!is.na(estimate), !is.na(se_estimate)) %>%
    dplyr::group_by(T, estimand, method_id, method, family) %>%
    dplyr::summarise(
      nsim = dplyr::n(),
      mc_sd = stats::sd(estimate, na.rm = TRUE),
      mean_se = mean(se_estimate, na.rm = TRUE),
      diff = mean_se - mc_sd,
      ratio = dplyr::if_else(mc_sd > 0, mean_se / mc_sd, NA_real_),
      .groups = "drop"
    )

# --------------------------------------------------------------------------------------------------
  # 7F. ML summaries
# --------------------------------------------------------------------------------------------------
  mse_df <- summarise_mean_mcse(
    data = ml_long %>% dplyr::filter(metric == "MSE"),
    value_col = "value",
    mean_name = "mean_mse",
    mcse_name = "mcse_mean_mse",
    group_cols = c("T", "target", "ml_method_id", "ml_method")
  )

  r2_df <- summarise_mean_mcse(
    data = ml_long %>% dplyr::filter(metric == "R2"),
    value_col = "value",
    mean_name = "mean_r2",
    mcse_name = "mcse_mean_r2",
    group_cols = c("T", "target", "ml_method_id", "ml_method")
  )

# --------------------------------------------------------------------------------------------------
  # 7G. Theoretical reference lines for MSE and R²
  #
  # If variance is 1, and true confounder-induced R² is known, then the best
  # achievable MSE is:
  #
  #   theoretical minimum MSE = 1 - true_R²
  #
  # For R², the reference line is the true confounder-induced R² itself.
# --------------------------------------------------------------------------------------------------
  theory_df <- df_all %>%
    dplyr::group_by(T) %>%
    dplyr::summarise(
      true_r2_x = mean(true_r2_x, na.rm = TRUE),
      true_r2_y = mean(true_r2_y, na.rm = TRUE),
      theoretical_min_mse_x = 1 - true_r2_x,
      theoretical_min_mse_y = 1 - true_r2_y,
      .groups = "drop"
    )

# --------------------------------------------------------------------------------------------------
  # 7H. Flag summaries
# --------------------------------------------------------------------------------------------------
  flag0_df <- run_level %>%
    dplyr::group_by(method_id, method, family) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      prop_flag0 = mean(flag0, na.rm = TRUE),
      .groups = "drop"
    )

  flag1_df <- run_level %>%
    dplyr::group_by(method_id, method, family) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      prop_flag1 = mean(flag1, na.rm = TRUE),
      .groups = "drop"
    )

  flag2_df <- run_level %>%
    dplyr::group_by(method_id, method, family) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      prop_flag2 = mean(flag2, na.rm = TRUE),
      .groups = "drop"
    )

# --------------------------------------------------------------------------------------------------
  # 7I. Improper-reason summary
# --------------------------------------------------------------------------------------------------
  improper_reason_df <- run_level %>%
    dplyr::filter(flag2 == 1 | !is.na(improper_reason)) %>%
    dplyr::mutate(
      improper_reason = dplyr::if_else(
        is.na(improper_reason),
        "unspecified improper reason",
        improper_reason
      )
    ) %>%
    dplyr::group_by(method_id, method, family, improper_reason) %>%
    dplyr::summarise(
      n_reason = dplyr::n(),
      .groups = "drop_last"
    ) %>%
    dplyr::mutate(
      n_improper_total = sum(n_reason),
      prop_reason = n_reason / n_improper_total
    ) %>%
    dplyr::ungroup()

  improper_reason_top_df <- improper_reason_df %>%
    dplyr::group_by(method_id, method, family) %>%
    dplyr::slice_max(order_by = n_reason, n = improper_top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

# ---- 8. ordering and palettes --------------------------------------------------------------------
  # Consistent ordering makes it easier to compare plots.

  method_key <- df_all %>%
    dplyr::distinct(method_id, method, family) %>%
    dplyr::arrange(family, method)

  method_levels <- method_key$method

  pal_method <- setNames(
    scales::hue_pal()(length(method_levels)),
    method_levels
  )

  ml_levels <- c("LM", "EN", "XGB")
  ml_levels_present <- intersect(ml_levels, unique(df_all$ml_method))
  pal_ml <- setNames(
    scales::hue_pal()(length(ml_levels_present)),
    ml_levels_present
  )

  relbias_df$method <- factor(as.character(relbias_df$method), levels = method_levels)
  rmse_df$method <- factor(as.character(rmse_df$method), levels = method_levels)
  se_df$method <- factor(as.character(se_df$method), levels = method_levels)
  detect_df$method <- factor(as.character(detect_df$method), levels = method_levels)
  se_check$method <- factor(as.character(se_check$method), levels = method_levels)

  flag0_df$method <- factor(as.character(flag0_df$method), levels = method_levels)
  flag1_df$method <- factor(as.character(flag1_df$method), levels = method_levels)
  flag2_df$method <- factor(as.character(flag2_df$method), levels = method_levels)
  improper_reason_df$method <- factor(as.character(improper_reason_df$method), levels = method_levels)
  improper_reason_top_df$method <- factor(as.character(improper_reason_top_df$method), levels = method_levels)

  mse_df$ml_method <- factor(as.character(mse_df$ml_method), levels = ml_levels_present)
  r2_df$ml_method <- factor(as.character(r2_df$ml_method), levels = ml_levels_present)

# ---- 9. family-specific sem overview plots -------------------------------------------------------
  # These are the main SEM overview plots.

  plot_relbias_clpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Relative bias",
    zero_line = TRUE
  )

  plot_relbias_riclpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "RI-CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Relative bias",
    zero_line = TRUE
  )

  plot_relbias_dpm <- make_family_plot(
    data = relbias_df,
    y = "rel_bias",
    y_se = "mcse_rel_bias",
    family_name = "DPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Relative bias",
    zero_line = TRUE
  )

  plot_rmse_clpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "RMSE"
  )

  plot_rmse_riclpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "RI-CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "RMSE"
  )

  plot_rmse_dpm <- make_family_plot(
    data = rmse_df,
    y = "rmse",
    y_se = "mcse_rmse",
    family_name = "DPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "RMSE"
  )

  plot_se_clpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Mean estimated SE"
  )

  plot_se_riclpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "RI-CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Mean estimated SE"
  )

  plot_se_dpm <- make_family_plot(
    data = se_df,
    y = "mean_se",
    y_se = "mcse_mean_se",
    family_name = "DPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Mean estimated SE"
  )

# ---- 10. family-specific power plots -------------------------------------------------------------
  # These are split by SEM family as requested.

  plot_power_clpm <- make_family_plot(
    data = detect_df,
    y = "detect_prob",
    y_se = "mcse_detect",
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
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
    palette = pal_method,
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
    palette = pal_method,
    y_label = paste0("Detection probability (alpha = ", alpha, ")"),
    ref_line = alpha,
    percent_y = TRUE,
    clamp_01 = TRUE
  )

# ---- 11. family-specific se calibration plots ----------------------------------------------------
  # Axis titles cleaned up as requested.

  plot_se_ratio_clpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE / Monte Carlo SD",
    ref_line = 1
  )

  plot_se_ratio_riclpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "RI-CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE / Monte Carlo SD",
    ref_line = 1
  )

  plot_se_ratio_dpm <- make_family_plot(
    data = se_check,
    y = "ratio",
    y_se = NULL,
    family_name = "DPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE / Monte Carlo SD",
    ref_line = 1
  )

  plot_se_diff_clpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE - Monte Carlo SD",
    ref_line = 0
  )

  plot_se_diff_riclpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "RI-CLPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE - Monte Carlo SD",
    ref_line = 0
  )

  plot_se_diff_dpm <- make_family_plot(
    data = se_check,
    y = "diff",
    y_se = NULL,
    family_name = "DPM",
    method_key = method_key,
    palette = pal_method,
    y_label = "Est. SE - Monte Carlo SD",
    ref_line = 0
  )

# ---- 12. ml overview plots -----------------------------------------------------------------------
  # Separate X and Y plots, each with the requested theoretical reference line.

  plot_mse_x <- make_line_plot(
    data = mse_df %>% dplyr::filter(target == "X"),
    y = "mean_mse",
    y_se = "mcse_mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF MSE",
    ref_df = theory_df,
    ref_y = "theoretical_min_mse_x"
  )

  plot_mse_y <- make_line_plot(
    data = mse_df %>% dplyr::filter(target == "Y"),
    y = "mean_mse",
    y_se = "mcse_mean_mse",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF MSE",
    ref_df = theory_df,
    ref_y = "theoretical_min_mse_y"
  )

  plot_r2_x <- make_line_plot(
    data = r2_df %>% dplyr::filter(target == "X"),
    y = "mean_r2",
    y_se = "mcse_mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF R²",
    ref_df = theory_df,
    ref_y = "true_r2_x"
  )

  plot_r2_y <- make_line_plot(
    data = r2_df %>% dplyr::filter(target == "Y"),
    y = "mean_r2",
    y_se = "mcse_mean_r2",
    color_var = "ml_method",
    group_var = "ml_method_id",
    palette = pal_ml,
    y_label = "Mean OOF R²",
    ref_df = theory_df,
    ref_y = "true_r2_y"
  )

# ---- 12b. true-delta plot ------------------------------------------------------------------------
  # This plot is based only on the stored data-generating deltas.
  #
  # Facets:
  #   one facet for each outcome ("x" and "y", if both are present)
  #
  # Lines:
  #   one line per confounder
  #
  # In constant scenarios this will show flat lines, which is often exactly what
  # you want to verify visually.

  if (nrow(true_delta_df) > 0) {
    plot_true_delta <- ggplot2::ggplot(
      true_delta_df,
      ggplot2::aes(
        x = T,
        y = delta,
        color = confounder,
        group = confounder
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~ outcome, nrow = 1, scales = "fixed") +
      ggplot2::scale_x_continuous(breaks = sort(unique(true_delta_df$T))) +
      ggplot2::labs(
        x = "Occasion",
        y = "True delta",
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 13),
        axis.text = ggplot2::element_text(size = 12),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 11)
      )
  } else {
    plot_true_delta <- NULL
  }

# ---- 13. flag bar plots --------------------------------------------------------------------------
  # Ordering is chosen for quick visual diagnosis.
  #
  # flag0:
  #   lowest proportion of proper runs on top
  #
  # flag1 / flag2:
  #   highest issue proportion on top

  plot_flag0 <- make_flag_bar_plot(
    data = flag0_df,
    x = "prop_flag0",
    x_label = "Proper runs",
    decreasing = FALSE,
    palette = pal_method
  )

  plot_flag1 <- make_flag_bar_plot(
    data = flag1_df,
    x = "prop_flag1",
    x_label = "Non-converged runs",
    decreasing = TRUE,
    palette = pal_method
  )

  plot_flag2 <- make_flag_bar_plot(
    data = flag2_df,
    x = "prop_flag2",
    x_label = "Improper runs",
    decreasing = TRUE,
    palette = pal_method
  )

# ---- 14. improper-reason plot --------------------------------------------------------------------
  # Keep the top improper reasons per method and facet by method.

  improper_reason_plot_df <- improper_reason_top_df %>%
    dplyr::mutate(
      reason_plot = reorder_within_base(
        x = improper_reason,
        by = prop_reason,
        within = method
      )
    )

  plot_improper_reasons <- ggplot2::ggplot(
    improper_reason_plot_df,
    ggplot2::aes(
      x = prop_reason,
      y = reason_plot,
      fill = method
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) strip_reorder_suffix(x)) +
    ggplot2::scale_x_continuous(labels = function(v) scales::percent(v, accuracy = 1)) +
    ggplot2::labs(
      x = "Share within improper runs",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::scale_fill_manual(values = pal_method, drop = FALSE)

# ---- 15. return ----------------------------------------------------------------------------------
  # Everything is returned as a named list so you can extract exactly what you
  # need without rerunning the function.

  list(
# --------------------------------------------------------------------------------------------------
    # Summary data frames
# --------------------------------------------------------------------------------------------------
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
    improper_reason_df = improper_reason_df,

# --------------------------------------------------------------------------------------------------
    # Family-specific SEM overview plots
# --------------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------------
    # Family-specific SE calibration plots
# --------------------------------------------------------------------------------------------------
    plot_se_ratio_clpm = plot_se_ratio_clpm,
    plot_se_ratio_riclpm = plot_se_ratio_riclpm,
    plot_se_ratio_dpm = plot_se_ratio_dpm,
    plot_se_diff_clpm = plot_se_diff_clpm,
    plot_se_diff_riclpm = plot_se_diff_riclpm,
    plot_se_diff_dpm = plot_se_diff_dpm,

# --------------------------------------------------------------------------------------------------
    # ML overview plots
# --------------------------------------------------------------------------------------------------
    plot_mse_x = plot_mse_x,
    plot_mse_y = plot_mse_y,
    plot_r2_x = plot_r2_x,
    plot_r2_y = plot_r2_y,

# --------------------------------------------------------------------------------------------------
    # Truth-only delta plot
# --------------------------------------------------------------------------------------------------
    plot_true_delta = plot_true_delta,

# --------------------------------------------------------------------------------------------------
    # Flag / diagnostic plots
# --------------------------------------------------------------------------------------------------
    plot_flag0 = plot_flag0,
    plot_flag1 = plot_flag1,
    plot_flag2 = plot_flag2,
    plot_improper_reasons = plot_improper_reasons
  )
}