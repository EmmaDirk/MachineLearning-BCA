# script to create the plots for simulation study 1
# -----------------------------------------------------------------------------

plot_sim_study_results <- function(results_sim, true_A) {

  #####################################################
  # true parameters from the data-generating model
  #####################################################

  true_params <- tibble(
    estimand = c("XY", "YX", "AX", "AY"),
    true     = c(
      true_A[2, 1],  # XY
      true_A[1, 2],  # YX
      true_A[1, 1],  # AX
      true_A[2, 2]   # AY
    )
  )

  #####################################################
  # build relative bias data frame
  #####################################################

  relbias_df <- results_sim %>%

    # only occasions with lagged parameters
    filter(occasion >= 2) %>%

    # reshape estimates to long format
    pivot_longer(
      cols = matches("^est(XY|YX|A|AY)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estXY_") ~ "XY",
        str_detect(param, "^estYX_") ~ "YX",
        str_detect(param, "^estA_")  ~ "AX",
        str_detect(param, "^estAY_") ~ "AY"
      ),

      method = param %>%
        str_remove("^est(XY|YX|A|AY)_") %>%
        recode(
          CLPM      = "CLPM",
          CLPM_Adj  = "True model",
          CLPM_LBCA = "CLPM linear BCA",
          DPM       = "DPM",
          RI_CLPM   = "RI-CLPM"
        )
    ) %>%

    # attach true parameter values
    left_join(true_params, by = "estimand") %>%

    # compute relative bias and Monte Carlo SE
    group_by(scenario, occasion, estimand, method) %>%
    reframe(
      nsim = n(),

      mean_est = mean(estimate),

      # relative bias
      rel_bias = (mean_est - true) / true,

      # Monte Carlo SE of relative bias
      mcse_rel_bias =
        sqrt(
          sum((estimate - mean_est)^2) /
            (nsim * (nsim - 1))
        ) / abs(true)
    ) %>%

    # set clean factor ordering
    mutate(
      scenario = factor(scenario, levels = c("constant", "stepwise")),
      method   = factor(
        method,
        levels = c(
          "CLPM",
          "True model",
          "CLPM linear BCA",
          "DPM",
          "RI-CLPM"
        )
      )
    )

  #####################################################
  # define colour palette for methods
  #####################################################

  pal_method <- viridis(
    n = nlevels(relbias_df$method),
    option = "viridis",
    begin  = 0.1,
    end    = 0.9
  )
  names(pal_method) <- levels(relbias_df$method)

  #####################################################
  # Plot 1: cross-lagged effect of X on Y (XY)
  #####################################################

  plot_relbias_XY <- relbias_df %>%
    filter(estimand == "XY") %>%

    ggplot(aes(
      x     = occasion,
      y     = rel_bias,
      color = method,
      group = method
    )) +

    # relative bias trajectories
    geom_line(linewidth = 0.9) +

    # add points at each occasion
    geom_point(size = 2) +

    # Monte Carlo SE error bars
    geom_errorbar(
      aes(
        ymin = rel_bias - mcse_rel_bias,
        ymax = rel_bias + mcse_rel_bias
      ),
      width = 0.15,
      linewidth = 0.4
    ) +

    # zero-bias reference line
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +

    # facet by scenario
    ggh4x::facet_wrap2(
      ~ scenario,
      nrow = 1,
      axes = "y"
    ) +

    # colours
    scale_color_manual(values = pal_method) +

    # labels
    labs(
      x = "Occasion",
      y = "Relative bias",
      color = NULL
    ) +

    # scientific theme
    theme_classic(base_size = 13) +
    theme(
      panel.spacing.x = unit(1.2, "cm"),
      strip.text      = element_text(size = 14),
      axis.title      = element_text(size = 13),
      axis.text       = element_text(size = 12),
      legend.position = "bottom",
      legend.text     = element_text(size = 12)
    )

  #####################################################
  # Build RMSE data for all lagged parameters
  #####################################################

  # true parameters from the data-generating model
  true_params <- tibble(
    estimand = c("XY", "YX", "AX", "AY"),
    true     = c(
      true_A[2, 1],  # XY
      true_A[1, 2],  # YX (true = 0)
      true_A[1, 1],  # AX
      true_A[2, 2]   # AY
    )
  )

  # build RMSE data frame
  rmse_df <- results_sim %>%

    # only occasions with lagged parameters
    filter(occasion >= 2) %>%

    # reshape estimates to long format
    pivot_longer(
      cols = matches("^est(XY|YX|A|AY)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estXY_") ~ "XY",
        str_detect(param, "^estYX_") ~ "YX",
        str_detect(param, "^estA_")  ~ "AX",
        str_detect(param, "^estAY_") ~ "AY"
      ),

      method = param %>%
        str_remove("^est(XY|YX|A|AY)_") %>%
        recode(
          CLPM      = "CLPM",
          CLPM_Adj  = "True model",
          CLPM_LBCA = "CLPM linear BCA",
          DPM       = "DPM",
          RI_CLPM   = "RI-CLPM"
        )
    ) %>%

    # attach true parameter values
    left_join(true_params, by = "estimand") %>%

    # compute RMSE and Monte Carlo SE
    group_by(scenario, occasion, estimand, method) %>%
    reframe(
      nsim = n(),

      # RMSE
      rmse = sqrt(mean((estimate - true)^2)),

      # Monte Carlo SE of RMSE
      mcse_rmse =
        sqrt(
          sum(
            ((estimate - true)^2 - mean((estimate - true)^2))^2
          ) /
            (nsim * (nsim - 1))
        ) / (2 * rmse)
    ) %>%

    # clean factor ordering
    mutate(
      scenario = factor(scenario, levels = c("constant", "stepwise")),
      method   = factor(
        method,
        levels = c(
          "CLPM",
          "True model",
          "CLPM linear BCA",
          "DPM",
          "RI-CLPM"
        )
      )
    )

  #####################################################
  # Colour palette for methods (shared across plots)
  #####################################################

  pal_method <- viridis(
    n = nlevels(rmse_df$method),
    option = "viridis",
    begin  = 0.1,
    end    = 0.9
  )
  names(pal_method) <- levels(rmse_df$method)

  #####################################################
  # Plot 1: cross-lagged effect of X on Y (XY)
  #####################################################

  plot_rmse_XY <- rmse_df %>%
    filter(estimand == "XY") %>%

    ggplot(aes(
      x     = occasion,
      y     = rmse,
      color = method,
      group = method
    )) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = rmse - mcse_rmse,
        ymax = rmse + mcse_rmse
      ),
      width = 0.15,
      linewidth = 0.4
    ) +
    ggh4x::facet_wrap2(~ scenario, nrow = 1, axes = "y") +
    scale_color_manual(values = pal_method) +
    labs(
      x = "Occasion",
      y = "RMSE",
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      panel.spacing.x = unit(1.2, "cm"),
      strip.text      = element_text(size = 14),
      axis.title      = element_text(size = 13),
      axis.text       = element_text(size = 12),
      legend.position = "bottom",
      legend.text     = element_text(size = 12)
    )

  #####################################################
  # create legendless plots
  #####################################################

  plot_relbias_XY_noleg <- plot_relbias_XY +
    theme(legend.position = "none")

  #####################################################
  # combine the two
  #####################################################

  combined_XY <- plot_relbias_XY_noleg / plot_rmse_XY

  return(list(
    combined_XY        = combined_XY,
    plot_relbias_XY    = plot_relbias_XY,
    plot_rmse_XY       = plot_rmse_XY,
    relbias_df         = relbias_df,
    rmse_df            = rmse_df,
    pal_method         = pal_method
  ))
}