# ============================================================
# SELF-CONTAINED QUICK DIAGNOSTIC FOR THE BCA SE PROBLEM
# ============================================================

library(here)

# ------------------------------------------------------------
# source your study scripts
# ------------------------------------------------------------
source(here("01_scripts", "02_development", "04_thesis_v1", "00_packages.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "01_delta_sampler.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "02_delta_trajectory.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "04_lavaan_model_string_builder.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "05_residualisers.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "06_model_fitters.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "07_bootstrap_helpers.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "08_fit_stat_extractors.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "09_one_replication_wrapper.R"))
source(here("01_scripts", "02_development", "04_thesis_v1", "10_simulation_function.R"))

# ============================================================
# SETTINGS
# ============================================================

set.seed(1233)

# ----- DGM -----
k <- 2
T_waves <- 5

Omega11 <- diag(2)

Delta1 <- sample_delta_1(
  k = k,
  Omega11 = Omega11,
  R2_total = 0.15,
  R2_interaction = 0,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1,
  force_positive = TRUE
)

Delta_list <- generate_Delta_constant(
  Delta1 = Delta1$Delta,
  T = T_waves
)

Phi <- matrix(c(
  0.20, 0.00,
  0.10, 0.20
), nrow = 2, byrow = TRUE)

Sigma <- matrix(c(
  1.00, 0.30,
  0.30, 1.00
), nrow = 2, byrow = TRUE)

# ----- model choices -----
diag_model <- list(
  N = 300,
  T = T_waves,
  k = k,
  residualizer = "linear",
  sem_model = "clpm",
  confounder_order = 1,
  exclude = NULL,
  free_loadings = FALSE,
  xgb_tuning = NULL
)

# ----- run settings -----
diag_run <- list(
  B = 60,
  reps = 8,
  occasion = 2,
  sim_seed_base = 5000,
  boot_seed_base = 9000
)

# residualizer args
residualizer_args_linear <- list(
  oof_folds = 5,
  seed = 123
)

residualizer_args_xgb <- list(
  oof_folds = 5,
  nthread = 1,
  seed = 123
)

residualizer_args <- if (diag_model$residualizer == "linear") {
  residualizer_args_linear
} else {
  residualizer_args_xgb
}

# ============================================================
# HELPERS
# ============================================================

get_sim_data <- function(sim_obj) {
  if (is.list(sim_obj) && !is.data.frame(sim_obj) && !is.null(sim_obj$data)) {
    sim_obj$data
  } else {
    sim_obj
  }
}

fit_once_extract <- function(
  df, T, k, residualizer, sem_model,
  confounder_order, exclude, free_loadings,
  xgb_tuning = NULL, residualizer_args = list()
) {
  fit_out <- fit_analysis_pipeline(
    df = df,
    T = T,
    k = k,
    residualizer = residualizer,
    sem_model = sem_model,
    confounder_order = confounder_order,
    exclude = exclude,
    free_loadings = free_loadings,
    xgb_tuning = xgb_tuning,
    residualizer_args = residualizer_args
  )

  lag <- extract_lagged_estimates(
    fit = fit_out$fit,
    T = T,
    model_type = sem_model
  )

  list(
    fit = fit_out$fit,
    lag = lag,
    flag = classify_fit_flag(fit_out$fit)
  )
}

collect_bootstrap_draws <- function(
  df, T, k, residualizer, sem_model,
  confounder_order, exclude, free_loadings,
  B, boot_seed, occasion,
  residualizer_args = list(),
  xgb_tuning = NULL,
  add_id_orig = TRUE,
  fresh_fold_seed = FALSE
) {

  if (!is.null(boot_seed)) set.seed(boot_seed)

  df_work <- as.data.frame(df)

  if (add_id_orig && !(".id_orig" %in% names(df_work))) {
    df_work$.id_orig <- seq_len(nrow(df_work))
  }

  out <- vector("list", B)

  for (b in seq_len(B)) {

    idx <- sample.int(
      n = nrow(df_work),
      size = nrow(df_work),
      replace = TRUE
    )

    df_b <- df_work[idx, , drop = FALSE]

    res_args_b <- residualizer_args

    if (fresh_fold_seed) {
      res_args_b$seed <- boot_seed + b
    }

    fit_b <- fit_analysis_pipeline(
      df = df_b,
      T = T,
      k = k,
      residualizer = residualizer,
      sem_model = sem_model,
      confounder_order = confounder_order,
      exclude = exclude,
      free_loadings = free_loadings,
      xgb_tuning = xgb_tuning,
      residualizer_args = res_args_b
    )

    lag_b <- extract_lagged_estimates(
      fit = fit_b$fit,
      T = T,
      model_type = sem_model
    )

    out[[b]] <- tibble::tibble(
      b = b,
      flag = classify_fit_flag(fit_b$fit),
      ARX = lag_b$ARX[occasion],
      ARY = lag_b$ARY[occasion],
      CXY = lag_b$CXY[occasion],
      CYX = lag_b$CYX[occasion]
    )
  }

  dplyr::bind_rows(out)
}

summarise_boot_draws <- function(draws_df, label) {
  draws_df %>%
    tidyr::pivot_longer(
      cols = c(ARX, ARY, CXY, CYX),
      names_to = "estimand",
      values_to = "estimate"
    ) %>%
    dplyr::group_by(estimand) %>%
    dplyr::summarise(
      variant = label,
      B_total = dplyr::n(),
      B_ok = sum(flag == 0 & !is.na(estimate)),
      prop_ok = mean(flag == 0 & !is.na(estimate)),
      boot_mean = mean(estimate[flag == 0], na.rm = TRUE),
      boot_sd = stats::sd(estimate[flag == 0], na.rm = TRUE),
      .groups = "drop"
    )
}

make_one_dataset <- function(seed) {
  sim <- simulate_panel_data(
    N = diag_model$N,
    T = diag_model$T,
    Phi = Phi,
    Delta_list = Delta_list,
    Omega11 = Omega11,
    Sigma = Sigma,
    seed = seed
  )
  get_sim_data(sim)
}

# ============================================================
# PART 1: ONE-DATASET CHECK
# ============================================================

df1 <- make_one_dataset(diag_run$sim_seed_base + 1)

df1_id <- df1
df1_id$.id_orig <- seq_len(nrow(df1_id))

main_current <- fit_once_extract(
  df = df1,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  xgb_tuning = diag_model$xgb_tuning,
  residualizer_args = residualizer_args
)

main_with_id <- fit_once_extract(
  df = df1_id,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  xgb_tuning = diag_model$xgb_tuning,
  residualizer_args = residualizer_args
)

main_compare <- tibble::tibble(
  estimand = c("ARX", "ARY", "CXY", "CYX"),
  main_current = c(
    main_current$lag$ARX[diag_run$occasion],
    main_current$lag$ARY[diag_run$occasion],
    main_current$lag$CXY[diag_run$occasion],
    main_current$lag$CYX[diag_run$occasion]
  ),
  main_with_id = c(
    main_with_id$lag$ARX[diag_run$occasion],
    main_with_id$lag$ARY[diag_run$occasion],
    main_with_id$lag$CXY[diag_run$occasion],
    main_with_id$lag$CYX[diag_run$occasion]
  )
) %>%
  dplyr::mutate(abs_diff = abs(main_current - main_with_id))

print(main_compare)

boot_current <- collect_bootstrap_draws(
  df = df1,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  B = diag_run$B,
  boot_seed = diag_run$boot_seed_base + 1,
  occasion = diag_run$occasion,
  residualizer_args = residualizer_args,
  xgb_tuning = diag_model$xgb_tuning,
  add_id_orig = TRUE,
  fresh_fold_seed = FALSE
)

boot_no_id <- collect_bootstrap_draws(
  df = df1,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  B = diag_run$B,
  boot_seed = diag_run$boot_seed_base + 1,
  occasion = diag_run$occasion,
  residualizer_args = residualizer_args,
  xgb_tuning = diag_model$xgb_tuning,
  add_id_orig = FALSE,
  fresh_fold_seed = FALSE
)

boot_fresh_seed <- collect_bootstrap_draws(
  df = df1,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  B = diag_run$B,
  boot_seed = diag_run$boot_seed_base + 1,
  occasion = diag_run$occasion,
  residualizer_args = residualizer_args,
  xgb_tuning = diag_model$xgb_tuning,
  add_id_orig = TRUE,
  fresh_fold_seed = TRUE
)

boot_no_id_fresh_seed <- collect_bootstrap_draws(
  df = df1,
  T = diag_model$T,
  k = diag_model$k,
  residualizer = diag_model$residualizer,
  sem_model = diag_model$sem_model,
  confounder_order = diag_model$confounder_order,
  exclude = diag_model$exclude,
  free_loadings = diag_model$free_loadings,
  B = diag_run$B,
  boot_seed = diag_run$boot_seed_base + 1,
  occasion = diag_run$occasion,
  residualizer_args = residualizer_args,
  xgb_tuning = diag_model$xgb_tuning,
  add_id_orig = FALSE,
  fresh_fold_seed = TRUE
)

one_dataset_summary <- dplyr::bind_rows(
  summarise_boot_draws(boot_current, "current"),
  summarise_boot_draws(boot_no_id, "no_id"),
  summarise_boot_draws(boot_fresh_seed, "fresh_fold_seed"),
  summarise_boot_draws(boot_no_id_fresh_seed, "no_id_plus_fresh_seed")
)

print(one_dataset_summary)

boot_plot_df <- dplyr::bind_rows(
  boot_current %>% dplyr::mutate(variant = "current"),
  boot_no_id %>% dplyr::mutate(variant = "no_id"),
  boot_fresh_seed %>% dplyr::mutate(variant = "fresh_fold_seed"),
  boot_no_id_fresh_seed %>% dplyr::mutate(variant = "no_id_plus_fresh_seed")
) %>%
  dplyr::filter(flag == 0) %>%
  tidyr::pivot_longer(
    cols = c(ARX, ARY, CXY, CYX),
    names_to = "estimand",
    values_to = "estimate"
  )

ggplot2::ggplot(boot_plot_df, ggplot2::aes(x = estimate)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(estimand ~ variant, scales = "free") +
  ggplot2::theme_classic(base_size = 12)

# ============================================================
# PART 2: MINI-SIMULATION CHECK
# ============================================================

run_mini_diag <- function(
  add_id_orig = TRUE,
  fresh_fold_seed = FALSE,
  label = "variant"
) {

  rep_out <- vector("list", diag_run$reps)

  for (r in seq_len(diag_run$reps)) {

    df_r <- make_one_dataset(diag_run$sim_seed_base + r)

    main_r <- fit_once_extract(
      df = df_r,
      T = diag_model$T,
      k = diag_model$k,
      residualizer = diag_model$residualizer,
      sem_model = diag_model$sem_model,
      confounder_order = diag_model$confounder_order,
      exclude = diag_model$exclude,
      free_loadings = diag_model$free_loadings,
      xgb_tuning = diag_model$xgb_tuning,
      residualizer_args = residualizer_args
    )

    boot_r <- collect_bootstrap_draws(
      df = df_r,
      T = diag_model$T,
      k = diag_model$k,
      residualizer = diag_model$residualizer,
      sem_model = diag_model$sem_model,
      confounder_order = diag_model$confounder_order,
      exclude = diag_model$exclude,
      free_loadings = diag_model$free_loadings,
      B = diag_run$B,
      boot_seed = diag_run$boot_seed_base + r,
      occasion = diag_run$occasion,
      residualizer_args = residualizer_args,
      xgb_tuning = diag_model$xgb_tuning,
      add_id_orig = add_id_orig,
      fresh_fold_seed = fresh_fold_seed
    )

    se_r <- summarise_boot_draws(boot_r, label) %>%
      dplyr::select(estimand, boot_sd)

    point_r <- tibble::tibble(
      R = r,
      estimand = c("ARX", "ARY", "CXY", "CYX"),
      estimate = c(
        main_r$lag$ARX[diag_run$occasion],
        main_r$lag$ARY[diag_run$occasion],
        main_r$lag$CXY[diag_run$occasion],
        main_r$lag$CYX[diag_run$occasion]
      )
    )

    rep_out[[r]] <- dplyr::left_join(point_r, se_r, by = "estimand")
  }

  dplyr::bind_rows(rep_out) %>%
    dplyr::group_by(estimand) %>%
    dplyr::summarise(
      variant = label,
      reps_used = sum(!is.na(estimate)),
      empirical_sd_across_reps = stats::sd(estimate, na.rm = TRUE),
      mean_boot_se = mean(boot_sd, na.rm = TRUE),
      ratio = mean_boot_se / empirical_sd_across_reps,
      .groups = "drop"
    )
}

mini_diag <- dplyr::bind_rows(
  run_mini_diag(TRUE,  FALSE, "current"),
  run_mini_diag(FALSE, FALSE, "no_id"),
  run_mini_diag(TRUE,  TRUE,  "fresh_fold_seed"),
  run_mini_diag(FALSE, TRUE,  "no_id_plus_fresh_seed")
)

print(mini_diag)

ggplot2::ggplot(mini_diag, ggplot2::aes(x = variant, y = ratio)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~ estimand, scales = "free_y") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Mean bootstrap SE / empirical SD across replications"
  )

# ============================================================
# OPTIONAL:
# diag_model$sem_model <- "riclpm"
# ============================================================