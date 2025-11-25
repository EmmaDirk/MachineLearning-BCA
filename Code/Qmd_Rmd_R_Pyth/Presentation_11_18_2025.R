# =========================================================
# Libraries
# =========================================================
library(tidyverse)   # dplyr, ggplot2, purrr, tibble, readr, etc.
library(ranger)      # random forest (for residuals)
library(lavaan)      # SEM (clpm / riclpm / dpm)
library(parallel)    # detectCores()
library(doParallel)  # registerDoParallel()
library(knitr)       # kable() for the table
library(xgboost)     # XGBoost
library(glmnet)      # LASSO

set.seed(123)

# =========================================================
# Helper: extract a specific effect from lavaan fit
# =========================================================
extract_effect <- function(fit, model_type, data_type, lhs, rhs) {
  pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)

  row <- pe %>%
    dplyr::filter(lhs == !!lhs, op == "~", rhs == !!rhs)

  if (nrow(row) == 0) {
    return(tibble::tibble(
      model_type = model_type,
      data_type  = data_type,
      lhs        = lhs,
      rhs        = rhs,
      est        = NA_real_,
      se         = NA_real_,
      z          = NA_real_,
      pvalue     = NA_real_,
      ci_lower   = NA_real_,
      ci_upper   = NA_real_,
      std_est    = NA_real_
    ))
  }

  row <- row[1, ]

  tibble::tibble(
    model_type = model_type,
    data_type  = data_type,
    lhs        = lhs,
    rhs        = rhs,
    est        = row$est,
    se         = row$se,
    z          = row$z,
    pvalue     = row$pvalue,
    ci_lower   = row$ci.lower,
    ci_upper   = row$ci.upper,
    std_est    = row$std.all
  )
}

# =========================================================
# Model specifications
# =========================================================

# -------------------------------------------------
# RI-CLPM model specification
# -------------------------------------------------
model_riclpm <- "
rix =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
riy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5
rix ~~ rix
riy ~~ riy
rix ~~ riy

x1 ~~ 0*x1; x2 ~~ 0*x2; x3 ~~ 0*x3; x4 ~~ 0*x4; x5 ~~ 0*x5
y1 ~~ 0*y1; y2 ~~ 0*y2; y3 ~~ 0*y3; y4 ~~ 0*y4; y5 ~~ 0*y5

wx1 =~ 1*x1; wx2 =~ 1*x2; wx3 =~ 1*x3; wx4 =~ 1*x4; wx5 =~ 1*x5
wy1 =~ 1*y1; wy2 =~ 1*y2; wy3 =~ 1*y3; wy4 =~ 1*y4; wy5 =~ 1*y5

rix ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
rix ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5
riy ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
riy ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5

wx1 ~~ wx1; wx2 ~~ wx2; wx3 ~~ wx3; wx4 ~~ wx4; wx5 ~~ wx5
wy1 ~~ wy1; wy2 ~~ wy2; wy3 ~~ wy3; wy4 ~~ wy4; wy5 ~~ wy5
wy1 ~~ wx1; wy2 ~~ wx2; wy3 ~~ wx3; wy4 ~~ wx4; wy5 ~~ wx5

wx2 ~ wx1 + wy1
wy2 ~ wx1 + wy1
wx3 ~ wx2 + wy2
wy3 ~ wx2 + wy2
wx4 ~ wx3 + wy3
wy4 ~ wx3 + wy3
wx5 ~ wx4 + wy4
wy5 ~ wx4 + wy4

x1 + x2 + x3 + x4 + x5 ~ mx*1
y1 + y2 + y3 + y4 + y5 ~ my*1
"

# -------------------------------------------------
# DPM model specification
# -------------------------------------------------
model_dpm <- "
fx =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
fy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5

fx ~~ x1 + y1
fy ~~ x1 + y1

x2 + y2 ~ x1 + y1
x3 + y3 ~ x2 + y2
x4 + y4 ~ x3 + y3
x5 + y5 ~ x4 + y4

x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
x4 ~~ y4
x5 ~~ y5

fx ~~ fx
fy ~~ fy
fx ~~ fy

x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
y1 ~~ y1
y2 ~~ y2
y3 ~~ y3
y4 ~~ y4
y5 ~~ y5

x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
y1 ~ 1
y2 ~ 1
y3 ~ 1
y4 ~ 1
y5 ~ 1
"

# -------------------------------------------------
# CLPM model specification
# -------------------------------------------------
model_clpm <- "
x2 + y2 ~ x1 + y1
x3 + y3 ~ x2 + y2
x4 + y4 ~ x3 + y3
x5 + y5 ~ x4 + y4

x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
x4 ~~ y4
x5 ~~ y5

x1 ~~ x1
y1 ~~ y1
x2 ~~ x2
y2 ~~ y2
x3 ~~ x3
y3 ~~ y3
x4 ~~ x4
y4 ~~ y4
x5 ~~ x5
y5 ~~ y5

x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
y1 ~ 1
y2 ~ 1
y3 ~ 1
y4 ~ 1
y5 ~ 1
"

# =========================================================
# Read data (6 datasets)
# =========================================================
linear      <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_linear_Pop.csv")
cubic       <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_cubic_Pop.csv")
exponential <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_exponential_Pop.csv")
sine        <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_Sin_Pop.csv")
tanh        <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_tanh_Pop.csv")
plateau     <- read_csv("C:/Users/Admin/Desktop/UU/Thesis umbrella/Thesis/Code/Data/Data_Sim_plateau_Pop.csv")

data_list <- list(
  linear      = linear,
  cubic       = cubic,
  exponential = exponential,
  sine        = sine,
  tanh        = tanh,
  plateau     = plateau
)

# =========================================================
# Small samples (n=2000) for residual plots: x1 ~ c1
# =========================================================
sample_small <- lapply(data_list, function(df) dplyr::sample_n(df, 2000))

sample_small <- lapply(sample_small, function(df) {
  mod <- lm(x1 ~ c1, data = df)
  df$resid <- resid(mod)
  df
})

purrr::imap(sample_small, ~ {
  ggplot(.x, aes(x = c1, y = resid)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    ggtitle(.y)
})

# =========================================================
# Big samples (n=10000)
# =========================================================
data_big <- lapply(data_list, function(df) dplyr::sample_n(df, 10000))

# =========================================================
# Residualizers: LM, tuned ranger, tuned XGBoost, LASSO
# =========================================================

# ---------------- LM residuals --------------------
add_lm_residuals <- function(df) {
  xs <- paste0("x", 1:5)
  for (x in xs) {
    f <- as.formula(paste0(x, " ~ c1 + c2 + c3"))
    mod <- lm(f, data = df)
    df[[paste0("res_lm_", x)]] <- resid(mod)
  }
  df
}

# -------------- Ranger residuals (tuned) ----------
add_ranger_residuals <- function(df) {
  xs <- paste0("x", 1:5)

  for (x in xs) {
    f <- as.formula(paste0(x, " ~ c1 + c2 + c3"))

    mtry_grid     <- 1:3
    min_node_grid <- c(1, 5, 10)
    best_err      <- Inf
    best_model    <- NULL

    for (mtry_val in mtry_grid) {
      for (min_node in min_node_grid) {
        mod_tmp <- ranger(
          formula        = f,
          data           = df,
          num.trees      = 200,
          mtry           = mtry_val,
          min.node.size  = min_node,
          importance     = "none",
          respect.unordered.factors = "order",
          num.threads    = parallel::detectCores() - 1
        )

        err <- mod_tmp$prediction.error

        if (!is.null(err) && err < best_err) {
          best_err   <- err
          best_model <- mod_tmp
        }
      }
    }

    preds <- predict(best_model, data = df)$predictions
    df[[paste0("res_rg_", x)]] <- df[[x]] - preds
  }

  df
}

# -------------- XGBoost residuals (tuned) ---------
add_xgb_residuals <- function(df) {
  xs <- paste0("x", 1:5)
  X  <- as.matrix(df[, c("c1", "c2", "c3")])

  for (x in xs) {
    y      <- df[[x]]
    dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

    cv <- xgboost::xgb.cv(
      data       = dtrain,
      nrounds    = 200,
      nfold      = 5,
      early_stopping_rounds = 10,
      objective  = "reg:squarederror",
      metrics    = "rmse",
      verbose    = 0,
      num.threads    = parallel::detectCores() - 1
    )

    best_nrounds <- cv$best_iteration

    mod <- xgboost::xgboost(
      data      = dtrain,
      objective = "reg:squarederror",
      nrounds   = best_nrounds,
      max_depth = 3,
      eta       = 0.1,
      verbose   = 0,
      num.threads    = parallel::detectCores() - 1
    )

    preds <- predict(mod, newdata = X)
    df[[paste0("res_xgb_", x)]] <- y - preds
  }

  df
}

# -------------- LASSO residuals (cv.glmnet) -------
add_lasso_residuals <- function(df) {
  xs <- paste0("x", 1:5)
  X  <- model.matrix(~ c1 + c2 + c3, data = df)[, -1, drop = FALSE]

  for (x in xs) {
    y <- df[[x]]

    fit <- glmnet::cv.glmnet(
      x      = X,
      y      = y,
      alpha  = 1,
      family = "gaussian",
      num.threads    = parallel::detectCores() - 1
    )

    preds <- predict(fit, newx = X, s = "lambda.min")[, 1]
    df[[paste0("res_lasso_", x)]] <- y - preds
  }

  df
}

# =========================================================
# Apply residualizers
# =========================================================
data_big_lm    <- lapply(data_big, add_lm_residuals)
data_big_rg    <- lapply(data_big, add_ranger_residuals)
data_big_xgb   <- lapply(data_big, add_xgb_residuals)
data_big_lasso <- lapply(data_big, add_lasso_residuals)

# =========================================================
# Build SEM-ready datasets for each residual type
# =========================================================
make_lm_residual_data <- function(df) {
  df %>%
    transmute(
      x1 = res_lm_x1,
      x2 = res_lm_x2,
      x3 = res_lm_x3,
      x4 = res_lm_x4,
      x5 = res_lm_x5,
      y1 = y1,
      y2 = y2,
      y3 = y3,
      y4 = y4,
      y5 = y5
    )
}

make_rg_residual_data <- function(df) {
  df %>%
    transmute(
      x1 = res_rg_x1,
      x2 = res_rg_x2,
      x3 = res_rg_x3,
      x4 = res_rg_x4,
      x5 = res_rg_x5,
      y1 = y1,
      y2 = y2,
      y3 = y3,
      y4 = y4,
      y5 = y5
    )
}

make_xgb_residual_data <- function(df) {
  df %>%
    transmute(
      x1 = res_xgb_x1,
      x2 = res_xgb_x2,
      x3 = res_xgb_x3,
      x4 = res_xgb_x4,
      x5 = res_xgb_x5,
      y1 = y1,
      y2 = y2,
      y3 = y3,
      y4 = y4,
      y5 = y5
    )
}

make_lasso_residual_data <- function(df) {
  df %>%
    transmute(
      x1 = res_lasso_x1,
      x2 = res_lasso_x2,
      x3 = res_lasso_x3,
      x4 = res_lasso_x4,
      x5 = res_lasso_x5,
      y1 = y1,
      y2 = y2,
      y3 = y3,
      y4 = y4,
      y5 = y5
    )
}

data_lm_sem    <- lapply(data_big_lm,    make_lm_residual_data)
data_rg_sem    <- lapply(data_big_rg,    make_rg_residual_data)
data_xgb_sem   <- lapply(data_big_xgb,   make_xgb_residual_data)
data_lasso_sem <- lapply(data_big_lasso, make_lasso_residual_data)

# =========================================================
# Fit CLPM, RI-CLPM, DPM on all datasets
# =========================================================
fit_all <- function(model_syntax, data_list) {
  lapply(data_list, function(df) {
    sem(
      model_syntax,
      data      = df,
      missing   = "fiml",
      estimator = "mlr"
    )
  })
}

# Raw data
fit_clpm_raw   <- fit_all(model_clpm,   data_big)
fit_riclpm_raw <- fit_all(model_riclpm, data_big)
fit_dpm_raw    <- fit_all(model_dpm,    data_big)

# CLPM on residualized data
fit_clpm_lm    <- fit_all(model_clpm, data_lm_sem)
fit_clpm_rg    <- fit_all(model_clpm, data_rg_sem)
fit_clpm_xgb   <- fit_all(model_clpm, data_xgb_sem)
fit_clpm_lasso <- fit_all(model_clpm, data_lasso_sem)

# =========================================================
# Collect fits into a single tibble
# =========================================================
collect_fits <- function(fit_list, model_type) {
  tibble(
    fit        = fit_list,
    model_type = model_type,
    data_type  = names(fit_list)
  )
}

fits_info <- bind_rows(
  collect_fits(fit_clpm_raw,   "clpm_raw"),
  collect_fits(fit_riclpm_raw, "riclpm_raw"),
  collect_fits(fit_dpm_raw,    "dpm_raw"),
  collect_fits(fit_clpm_lm,    "clpm_lm_resid"),
  collect_fits(fit_clpm_rg,    "clpm_rg_resid"),
  collect_fits(fit_clpm_xgb,   "clpm_xgb_resid"),
  collect_fits(fit_clpm_lasso, "clpm_lasso_resid")
) %>%
  mutate(
    lhs = dplyr::case_when(
      model_type == "riclpm_raw" ~ "wy2",
      TRUE                       ~ "y2"
    ),
    rhs = dplyr::case_when(
      model_type == "riclpm_raw" ~ "wx1",
      TRUE                       ~ "x1"
    )
  )

# =========================================================
# Extract x1 -> y2 effects from all models
# =========================================================
effects_x1_y2 <- purrr::pmap_dfr(
  .l = list(
    fit        = fits_info$fit,
    model_type = fits_info$model_type,
    data_type  = fits_info$data_type,
    lhs        = fits_info$lhs,
    rhs        = fits_info$rhs
  ),
  .f = extract_effect
)

effects_x1_y2_tbl <- effects_x1_y2 %>%
  arrange(model_type, data_type) %>%
  select(
    model_type, data_type, lhs, rhs,
    est, se, z, pvalue, ci_lower, ci_upper, std_est
  )

effects_x1_y2_tbl

effects_x1_y2_tbl %>%
  knitr::kable(
    digits  = 3,
    caption = "Estimated cross-lagged effect x1 → y2 across models and datasets"
  )

# =========================================================
# Bias calculation and plot
# =========================================================
true_value <- 0.1

bias_tbl <- effects_x1_y2_tbl %>%
  mutate(
    bias       = est - true_value,
    model_type = factor(model_type),
    data_type  = factor(data_type)
  )

bias_summary <- bias_tbl %>%
  group_by(model_type, data_type) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias   = sd(bias, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(
  bias_tbl,
  aes(
    x     = model_type,
    y     = bias,
    fill  = data_type,
    color = data_type
  )
) +
  geom_boxplot(
    position      = position_dodge(width = 0.8),
    alpha         = 0.4,
    outlier.shape = NA
  ) +
  geom_errorbar(
    data        = bias_summary,
    inherit.aes = FALSE,
    aes(
      x    = model_type,
      ymin = mean_bias - sd_bias,
      ymax = mean_bias + sd_bias,
      color = data_type
    ),
    width     = 0.15,
    position  = position_dodge(width = 0.8),
    linewidth = 0.6
  ) +
  geom_point(
    data        = bias_summary,
    inherit.aes = FALSE,
    aes(
      x     = model_type,
      y     = mean_bias,
      color = data_type
    ),
    size     = 2.8,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Bias of x1 → y2 across model and data types",
    x     = "Model type",
    y     = "Bias (estimate − 0.1)",
    color = "Data type",
    fill  = "Data type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

