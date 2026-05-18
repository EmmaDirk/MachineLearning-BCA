# This function is based on the plot_delta_trajectories function in
# 01_scripts/03_exploration/02_sim_engine/02_delta_trajectory/v1_testing_script.R
# and it converts a delta trajectory into a nice plot
# 
#  dependencies:
# library(ggplot2)
# library(viridis)

plot_delta_trajectories <- function(
    D_const,
    D_step,
    D_mix,
    scenario_names = c("constant", "stepwise", "stepwise mixture"),
    free_y = TRUE
) {

  # ------------------------------------------------------------
  # helper: convert one D_list into long format
  # ------------------------------------------------------------
  D_list_to_long <- function(D_list, scenario_name) {

    out <- vector("list", length(D_list))

    for (t in seq_along(D_list)) {
      D_t <- D_list[[t]]

      out[[t]] <- data.frame(
        time = t,
        outcome = rep(rownames(D_t), each = ncol(D_t)),
        coefficient = rep(colnames(D_t), times = nrow(D_t)),
        delta = as.vector(t(D_t)),
        scenario = scenario_name,
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, out)
  }

  # ------------------------------------------------------------
  # combine data
  # ------------------------------------------------------------
  plot_dat <- rbind(
    D_list_to_long(D_const, scenario_names[1]),
    D_list_to_long(D_step,  scenario_names[2]),
    D_list_to_long(D_mix,   scenario_names[3])
  )

  # explicit scenario order
  plot_dat$scenario <- factor(
    plot_dat$scenario,
    levels = scenario_names
  )

  # stable coefficient colors
  plot_dat$coefficient <- factor(
    plot_dat$coefficient,
    levels = unique(plot_dat$coefficient)
  )

  # facet scaling
  facet_scales <- if (free_y) "free_y" else "fixed"

  # ------------------------------------------------------------
  # build plot
  # ------------------------------------------------------------
  p <- ggplot(
    plot_dat,
    aes(x = time, y = delta, group = coefficient, color = coefficient)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_grid(outcome ~ scenario, scales = facet_scales) +
    labs(
      x = "Time",
      y = expression(delta),
      color = "Coefficient",
      title = "Delta trajectories under three scenarios"
    ) +
    scale_color_viridis_d(option = "D") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )

  return(p)
}