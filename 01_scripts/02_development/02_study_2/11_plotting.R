# script to create the plots for simulation study 2
#
# This function takes the simulation results data frame and the true A matrix from the data-generating model,
# and produces summary plots for the cross-lagged effect GammaXY (X <- Y), along with the underlying summary tables.
#
# 1) It defines the true parameters from the generating model:
# - GammaXY and GammaYX are taken from the off-diagonal entries of true_A
# - BetaX and BetaY are taken from the diagonal entries of true_A
#
# 2) It defines a recoding map to turn raw column suffixes (CLPM, RI_CLPM, etc.) into human-readable method labels.
#
# 3) It builds a relative-bias summary data frame:
# - it keeps occasions >= 2 (lagged parameters exist only from wave 2 onward)
# - it reshapes the estimate columns to long format
# - it identifies which estimand each estimate belongs to (GammaXY, GammaYX, BetaX, BetaY)
# - it computes relative bias and a Monte Carlo SE for the relative bias
#
# 4) It builds an RMSE summary data frame with the same reshaping and grouping logic:
# - it computes RMSE and a Monte Carlo SE for RMSE
#
# 5) It creates two plots for GammaXY:
# - relative bias over occasions (with Monte Carlo SE error bars)
# - RMSE over occasions (with Monte Carlo SE error bars)
#
# 6) It returns the plots, palettes, and the summary data frames for further use.
# -----------------------------------------------------------------------------

plot_sim_study_results <- function(results_sim, true_A) {

  #####################################################
  # true parameters from the data-generating model
  #####################################################

  true_params <- tibble(
    estimand = c("GammaXY", "GammaYX", "BetaX", "BetaY"),
    true     = c(
      true_A[2, 1],  # GammaXY
      true_A[1, 2],  # GammaYX
      true_A[1, 1],  # BetaX
      true_A[2, 2]   # BetaY
    )
  )

  #####################################################
  # method recoding map
  #####################################################

  method_map <- c(
    CLPM      = "CLPM",
    CLPM_Adj  = "True model",
    CLPM_LBCA = "CLPM linear BCA",
    CLPM_XGB  = "CLPM XGB",
    DPM       = "DPM",
    RI_CLPM   = "RI-CLPM"
  )

  #####################################################
  # build relative bias data frame
  #####################################################

  relbias_df <- results_sim %>%

    # only occasions with lagged parameters
    filter(occasion >= 2) %>%

    # reshape estimates to long format
    pivot_longer(
      cols = matches("^est(GammaXY|GammaYX|BetaX|BetaY)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estGammaXY_") ~ "GammaXY",
        str_detect(param, "^estGammaYX_") ~ "GammaYX",
        str_detect(param, "^estBetaX_")   ~ "BetaX",
        str_detect(param, "^estBetaY_")   ~ "BetaY"
      ),

      method_raw = param %>%
        str_remove("^est(GammaXY|GammaYX|BetaX|BetaY)_"),

      method = method_raw %>%
        recode(!!!method_map)
    ) %>%

    # keep only methods that exist in the results
    filter(!is.na(method)) %>%

    # attach true parameter values
    left_join(true_params, by = "estimand") %>%

    # compute relative bias and Monte Carlo SE
    group_by(scenario, occasion, estimand, method) %>%
    reframe(
      nsim = n(),

      mean_est = mean(estimate),

      rel_bias = (mean_est - true) / true,

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
          "CLPM XGB",
          "DPM",
          "RI-CLPM"
        )
      ),
      estimand = factor(estimand, levels = c("GammaXY", "GammaYX", "BetaX", "BetaY"))
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
  # Plot 1: relative bias for cross-lagged effect GammaXY
  #####################################################

  plot_relbias_GammaXY <- relbias_df %>%
    filter(estimand == "GammaXY") %>%

    ggplot(aes(
      x     = occasion,
      y     = rel_bias,
      color = method,
      group = method
    )) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = rel_bias - mcse_rel_bias,
        ymax = rel_bias + mcse_rel_bias
      ),
      width = 0.15,
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    ggh4x::facet_wrap2(
      ~ scenario,
      nrow = 1,
      axes = "y"
    ) +
    scale_color_manual(values = pal_method) +
    labs(
      x = "Occasion",
      y = "Relative bias",
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
  # build RMSE data frame
  #####################################################

  rmse_df <- results_sim %>%

    # only occasions with lagged parameters
    filter(occasion >= 2) %>%

    # reshape estimates to long format
    pivot_longer(
      cols = matches("^est(GammaXY|GammaYX|BetaX|BetaY)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estGammaXY_") ~ "GammaXY",
        str_detect(param, "^estGammaYX_") ~ "GammaYX",
        str_detect(param, "^estBetaX_")   ~ "BetaX",
        str_detect(param, "^estBetaY_")   ~ "BetaY"
      ),

      method_raw = param %>%
        str_remove("^est(GammaXY|GammaYX|BetaX|BetaY)_"),

      method = method_raw %>%
        recode(!!!method_map)
    ) %>%

    # keep only methods that exist in the results
    filter(!is.na(method)) %>%

    # attach true parameter values
    left_join(true_params, by = "estimand") %>%

    # compute RMSE and Monte Carlo SE
    group_by(scenario, occasion, estimand, method) %>%
    reframe(
      nsim = n(),

      rmse = sqrt(mean((estimate - true)^2)),

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
          "CLPM XGB",
          "DPM",
          "RI-CLPM"
        )
      ),
      estimand = factor(estimand, levels = c("GammaXY", "GammaYX", "BetaX", "BetaY"))
    )

  #####################################################
  # colour palette for methods
  #####################################################

  pal_method <- viridis(
    n = nlevels(rmse_df$method),
    option = "viridis",
    begin  = 0.1,
    end    = 0.9
  )
  names(pal_method) <- levels(rmse_df$method)

  #####################################################
  # Plot 2: RMSE for cross-lagged effect GammaXY
  #####################################################

  plot_rmse_GammaXY <- rmse_df %>%
    filter(estimand == "GammaXY") %>%

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

  plot_relbias_GammaXY_noleg <- plot_relbias_GammaXY +
    theme(legend.position = "none")

  #####################################################
  # combine the two
  #####################################################

  combined_GammaXY <- plot_relbias_GammaXY_noleg / plot_rmse_GammaXY

  return(list(
    combined_GammaXY        = combined_GammaXY,
    plot_relbias_GammaXY    = plot_relbias_GammaXY,
    plot_rmse_GammaXY       = plot_rmse_GammaXY,
    relbias_df              = relbias_df,
    rmse_df                 = rmse_df,
    pal_method              = pal_method
  ))
}
