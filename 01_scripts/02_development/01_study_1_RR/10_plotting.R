# This function takes the combined simulation output (results_sim) and the true 2x2 data-generating
# lag matrix (true_A), then produces summary performance plots for the cross-lagged X->Y effect:
# - Relative bias over time (with Monte Carlo SE error bars)
# - RMSE over time (with Monte Carlo SE error bars)
#
# Naming logic used throughout:
# - autoregressive effects are denoted beta:
#     beta_X = A[1,1] (X_t-1 -> X_t)
#     beta_Y = A[2,2] (Y_t-1 -> Y_t)
# - cross-lagged effects are denoted gamma:
#     gamma_XY = A[2,1] (X_t-1 -> Y_t)
#     gamma_YX = A[1,2] (Y_t-1 -> X_t)
#
# What the function does:
# 1) Builds a small "true_params" lookup table mapping each estimand (gamma_XY, gamma_YX, beta_X, beta_Y) to its true value.
# 2) Reshapes the simulation output from wide to long format for all lagged estimates (occasion >= 2).
# 3) Parses column names to identify which estimand (beta/gamma) and which method (CLPM, RI-CLPM, DPM, etc.) produced each estimate.
# 4) Computes, per scenario × occasion × estimand × method:
#    - relative bias and its Monte Carlo SE
#    - RMSE and its Monte Carlo SE
# 5) Creates two ggplots for gamma_XY (X -> Y): one for relative bias and one for RMSE, faceted by scenario.
# 6) Returns the plots plus the underlying summary data frames and the method colour palette.
# -----------------------------------------------------------------------------

plot_sim_study_results <- function(results_sim, true_A) {

  #####################################################
  # true parameters from the data-generating model
  #####################################################

  true_params <- tibble(
    estimand = c("gamma_XY", "gamma_YX", "beta_X", "beta_Y"),
    true     = c(
      true_A[2, 1],  # gamma_XY: X -> Y
      true_A[1, 2],  # gamma_YX: Y -> X
      true_A[1, 1],  # beta_X: X -> X
      true_A[2, 2]   # beta_Y: Y -> Y
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
      cols = matches("^est(BetaX|BetaY|GammaXY|GammaYX)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estGammaXY_") ~ "gamma_XY",
        str_detect(param, "^estGammaYX_") ~ "gamma_YX",
        str_detect(param, "^estBetaX_")   ~ "beta_X",
        str_detect(param, "^estBetaY_")   ~ "beta_Y"
      ),

      method = param %>%
        str_remove("^est(BetaX|BetaY|GammaXY|GammaYX)_") %>%
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
  # Plot 1: cross-lagged effect of X on Y (gamma_XY)
  #####################################################

  plot_relbias_gamma_XY <- relbias_df %>%
    filter(estimand == "gamma_XY") %>%

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

  # (reuse true_params defined above)

  rmse_df <- results_sim %>%

    # only occasions with lagged parameters
    filter(occasion >= 2) %>%

    # reshape estimates to long format
    pivot_longer(
      cols = matches("^est(BetaX|BetaY|GammaXY|GammaYX)_"),
      names_to  = "param",
      values_to = "estimate"
    ) %>%

    # drop failed fits
    filter(!is.na(estimate)) %>%

    # identify estimand type and method
    mutate(
      estimand = case_when(
        str_detect(param, "^estGammaXY_") ~ "gamma_XY",
        str_detect(param, "^estGammaYX_") ~ "gamma_YX",
        str_detect(param, "^estBetaX_")   ~ "beta_X",
        str_detect(param, "^estBetaY_")   ~ "beta_Y"
      ),

      method = param %>%
        str_remove("^est(BetaX|BetaY|GammaXY|GammaYX)_") %>%
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
  # Plot 2: RMSE for cross-lagged effect of X on Y (gamma_XY)
  #####################################################

  plot_rmse_gamma_XY <- rmse_df %>%
    filter(estimand == "gamma_XY") %>%

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

  plot_relbias_gamma_XY_noleg <- plot_relbias_gamma_XY +
    theme(legend.position = "none")

  #####################################################
  # combine the two
  #####################################################

  combined_gamma_XY <- plot_relbias_gamma_XY_noleg / plot_rmse_gamma_XY

  return(list(
    combined_gamma_XY        = combined_gamma_XY,
    plot_relbias_gamma_XY    = plot_relbias_gamma_XY,
    plot_rmse_gamma_XY       = plot_rmse_gamma_XY,
    relbias_df               = relbias_df,
    rmse_df                  = rmse_df,
    pal_method               = pal_method
  ))
}
