library(tidyverse)
library(ranger)
library(lavaan)

# parallel processing
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# reading data
linear    <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_linear_Pop.csv")
cubic_10  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.10_Pop.csv")
cubic_18  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.18_Pop.csv")
cubic_41  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.41_Pop.csv")
cubic_89  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.89_Pop.csv")

# take samples of each dataset
slinear <- linear %>% sample_n(2000)
scubic_10 <- cubic_10 %>% sample_n(2000)
scubic_18 <- cubic_18 %>% sample_n(2000)
scubic_41 <- cubic_41 %>% sample_n(2000)
scubic_89 <- cubic_89 %>% sample_n(2000)

# regress x1 on c1
mod_lin     <- lm(x1 ~ c1, data = slinear)
mod_cub_10  <- lm(x1 ~ c1, data = scubic_10)
mod_cub_18  <- lm(x1 ~ c1, data = scubic_18)
mod_cub_41  <- lm(x1 ~ c1, data = scubic_41)
mod_cub_89  <- lm(x1 ~ c1, data = scubic_89)

# add residuals as columns in the data
slinear$resid    <- resid(mod_lin)
scubic_10$resid  <- resid(mod_cub_10)
scubic_18$resid  <- resid(mod_cub_18)
scubic_41$resid  <- resid(mod_cub_41)
scubic_89$resid  <- resid(mod_cub_89)

# plot residuals versus predictor
ggplot(slinear, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Linear")

ggplot(scubic_10, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.10")

ggplot(scubic_18, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.18")

ggplot(scubic_41, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.41")

ggplot(scubic_89, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.89")

# take samples of each dataset
linear <- linear %>% sample_n(10000)
cubic_10 <- cubic_10 %>% sample_n(10000)
cubic_18 <- cubic_18 %>% sample_n(10000)
cubic_41 <- cubic_41 %>% sample_n(10000)
cubic_89 <- cubic_89 %>% sample_n(10000)
# linear
lm_lin_x1 <- lm(x1 ~ c1 + c2 + c3, data = linear)
lm_lin_x2 <- lm(x2 ~ c1 + c2 + c3, data = linear)
lm_lin_x3 <- lm(x3 ~ c1 + c2 + c3, data = linear)
lm_lin_x4 <- lm(x4 ~ c1 + c2 + c3, data = linear)
lm_lin_x5 <- lm(x5 ~ c1 + c2 + c3, data = linear)

# cubic_10
lm_cub10_x1 <- lm(x1 ~ c1 + c2 + c3, data = cubic_10)
lm_cub10_x2 <- lm(x2 ~ c1 + c2 + c3, data = cubic_10)
lm_cub10_x3 <- lm(x3 ~ c1 + c2 + c3, data = cubic_10)
lm_cub10_x4 <- lm(x4 ~ c1 + c2 + c3, data = cubic_10)
lm_cub10_x5 <- lm(x5 ~ c1 + c2 + c3, data = cubic_10)

# cubic_18
lm_cub18_x1 <- lm(x1 ~ c1 + c2 + c3, data = cubic_18)
lm_cub18_x2 <- lm(x2 ~ c1 + c2 + c3, data = cubic_18)
lm_cub18_x3 <- lm(x3 ~ c1 + c2 + c3, data = cubic_18)
lm_cub18_x4 <- lm(x4 ~ c1 + c2 + c3, data = cubic_18)
lm_cub18_x5 <- lm(x5 ~ c1 + c2 + c3, data = cubic_18)

# cubic_41
lm_cub41_x1 <- lm(x1 ~ c1 + c2 + c3, data = cubic_41)
lm_cub41_x2 <- lm(x2 ~ c1 + c2 + c3, data = cubic_41)
lm_cub41_x3 <- lm(x3 ~ c1 + c2 + c3, data = cubic_41)
lm_cub41_x4 <- lm(x4 ~ c1 + c2 + c3, data = cubic_41)
lm_cub41_x5 <- lm(x5 ~ c1 + c2 + c3, data = cubic_41)

# cubic_89
lm_cub89_x1 <- lm(x1 ~ c1 + c2 + c3, data = cubic_89)
lm_cub89_x2 <- lm(x2 ~ c1 + c2 + c3, data = cubic_89)
lm_cub89_x3 <- lm(x3 ~ c1 + c2 + c3, data = cubic_89)
lm_cub89_x4 <- lm(x4 ~ c1 + c2 + c3, data = cubic_89)
lm_cub89_x5 <- lm(x5 ~ c1 + c2 + c3, data = cubic_89)

# linear
linear <- linear %>%
  mutate(
    res_lm_x1 = resid(lm_lin_x1),
    res_lm_x2 = resid(lm_lin_x2),
    res_lm_x3 = resid(lm_lin_x3),
    res_lm_x4 = resid(lm_lin_x4),
    res_lm_x5 = resid(lm_lin_x5)
  )

# cubic_10
cubic_10 <- cubic_10 %>%
  mutate(
    res_lm_x1 = resid(lm_cub10_x1),
    res_lm_x2 = resid(lm_cub10_x2),
    res_lm_x3 = resid(lm_cub10_x3),
    res_lm_x4 = resid(lm_cub10_x4),
    res_lm_x5 = resid(lm_cub10_x5)
  )

# cubic_18
cubic_18 <- cubic_18 %>%
  mutate(
    res_lm_x1 = resid(lm_cub18_x1),
    res_lm_x2 = resid(lm_cub18_x2),
    res_lm_x3 = resid(lm_cub18_x3),
    res_lm_x4 = resid(lm_cub18_x4),
    res_lm_x5 = resid(lm_cub18_x5)
  )

# cubic_41
cubic_41 <- cubic_41 %>%
  mutate(
    res_lm_x1 = resid(lm_cub41_x1),
    res_lm_x2 = resid(lm_cub41_x2),
    res_lm_x3 = resid(lm_cub41_x3),
    res_lm_x4 = resid(lm_cub41_x4),
    res_lm_x5 = resid(lm_cub41_x5)
  )

# cubic_89
cubic_89 <- cubic_89 %>%
  mutate(
    res_lm_x1 = resid(lm_cub89_x1),
    res_lm_x2 = resid(lm_cub89_x2),
    res_lm_x3 = resid(lm_cub89_x3),
    res_lm_x4 = resid(lm_cub89_x4),
    res_lm_x5 = resid(lm_cub89_x5)
  )

# function to add ranger residuals
add_ranger_residuals <- function(df) {
  xs <- paste0("x", 1:5)

  for (x in xs) {
    f <- as.formula(paste0(x, " ~ c1 + c2 + c3"))
    mod <- ranger(f, data = df, num.trees = 100)  # reduce trees for speed if you like

    preds <- predict(mod, data = df)$predictions
    df[[paste0("res_rg_", x)]] <- df[[x]] - preds
  }

  df
}

linear_rg    <- add_ranger_residuals(linear)
cubic_10_rg  <- add_ranger_residuals(cubic_10)
cubic_18_rg  <- add_ranger_residuals(cubic_18)
cubic_41_rg  <- add_ranger_residuals(cubic_41)
cubic_89_rg  <- add_ranger_residuals(cubic_89)

# ri-clpm model specification
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

# dpm model specification
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

# clpm model specification
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

# fit a "normal" CLPM model to all data
fit_clpm_linear <- sem(model_clpm,
                       data = linear,
                       missing = "fiml",
                       estimator = "mlr")

fit_clpm_cubic_10 <- sem(model_clpm,
                         data = cubic_10,
                         missing = "fiml",
                         estimator = "mlr")

fit_clpm_cubic_18 <- sem(model_clpm,
                         data = cubic_18,
                         missing = "fiml",
                         estimator = "mlr")

fit_clpm_cubic_41 <- sem(model_clpm,
                         data = cubic_41,
                         missing = "fiml",
                         estimator = "mlr")

fit_clpm_cubic_89 <- sem(model_clpm,
                         data = cubic_89,
                         missing = "fiml",
                         estimator = "mlr")

# fit an RI-CLPM model to all data
fit_riclpm_linear <- sem(model_riclpm,
                         data = linear,
                         missing = "fiml",
                         estimator = "mlr")

fit_riclpm_cubic_10 <- sem(model_riclpm,
                           data = cubic_10,
                           missing = "fiml",
                           estimator = "mlr")

fit_riclpm_cubic_18 <- sem(model_riclpm,
                           data = cubic_18,
                           missing = "fiml",
                           estimator = "mlr")

fit_riclpm_cubic_41 <- sem(model_riclpm,
                           data = cubic_41,
                           missing = "fiml",
                           estimator = "mlr")

fit_riclpm_cubic_89 <- sem(model_riclpm,
                           data = cubic_89,
                           missing = "fiml",
                           estimator = "mlr")

# fit a DPM model to all data
fit_dpm_linear <- sem(model_dpm,
                      data = linear,
                      missing = "fiml",
                      estimator = "mlr")

fit_dpm_cubic_10 <- sem(model_dpm,
                        data = cubic_10,
                        missing = "fiml",
                        estimator = "mlr")

fit_dpm_cubic_18 <- sem(model_dpm,
                        data = cubic_18,
                        missing = "fiml",
                        estimator = "mlr")

fit_dpm_cubic_41 <- sem(model_dpm,
                        data = cubic_41,
                        missing = "fiml",
                        estimator = "mlr")

fit_dpm_cubic_89 <- sem(model_dpm,
                        data = cubic_89,
                        missing = "fiml",
                        estimator = "mlr")

# Fit a CLPM to all linear residual data
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

linear_lm_dat   <- make_lm_residual_data(linear)
cubic10_lm_dat  <- make_lm_residual_data(cubic_10)
cubic18_lm_dat  <- make_lm_residual_data(cubic_18)
cubic41_lm_dat  <- make_lm_residual_data(cubic_41)
cubic89_lm_dat  <- make_lm_residual_data(cubic_89)

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

linear_rg_dat   <- make_rg_residual_data(linear_rg)
cubic10_rg_dat  <- make_rg_residual_data(cubic_10_rg)
cubic18_rg_dat  <- make_rg_residual_data(cubic_18_rg)
cubic41_rg_dat  <- make_rg_residual_data(cubic_41_rg)
cubic89_rg_dat  <- make_rg_residual_data(cubic_89_rg)

# fit CLPM to lm residual data
fit_clpm_lm_linear <- sem(model_clpm,
                           data = linear_lm_dat,
                           missing = "fiml",
                           estimator = "mlr")

fit_clpm_lm_cubic_10 <- sem(model_clpm,
                             data = cubic10_lm_dat,
                             missing = "fiml",
                             estimator = "mlr")

fit_clpm_lm_cubic_18 <- sem(model_clpm,
                             data = cubic18_lm_dat,
                             missing = "fiml",
                             estimator = "mlr")

fit_clpm_lm_cubic_41 <- sem(model_clpm,
                             data = cubic41_lm_dat,
                             missing = "fiml",
                             estimator = "mlr")

fit_clpm_lm_cubic_89 <- sem(model_clpm,
                             data = cubic89_lm_dat,
                             missing = "fiml",
                             estimator = "mlr")

# fit CLPM to ranger residual data
fit_clpm_rg_linear <- sem(model_clpm,
                          data = linear_rg_dat,
                          missing = "fiml",
                          estimator = "mlr")

fit_clpm_rg_cubic_10 <- sem(model_clpm,
                            data = cubic10_rg_dat,
                            missing = "fiml",
                            estimator = "mlr")

fit_clpm_rg_cubic_18 <- sem(model_clpm,
                            data = cubic18_rg_dat,
                            missing = "fiml",
                            estimator = "mlr")

fit_clpm_rg_cubic_41 <- sem(model_clpm,
                            data = cubic41_rg_dat,
                            missing = "fiml",
                            estimator = "mlr")

fit_clpm_rg_cubic_89 <- sem(model_clpm,
                            data = cubic89_rg_dat,
                            missing = "fiml",
                            estimator = "mlr")

extract_effect <- function(fit, model_type, data_type, lhs, rhs) {
  pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  
  row <- pe %>% 
    filter(lhs == !!lhs, op == "~", rhs == !!rhs)
  
  if (nrow(row) == 0) {
    return(tibble(
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
  tibble(
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

# ---------------------------------------------
# list all models + which parameter to extract
# ---------------------------------------------
fits_info <- tibble::tibble(
  fit = list(
    # CLPM on raw data
    fit_clpm_linear,
    fit_clpm_cubic_10,
    fit_clpm_cubic_18,
    fit_clpm_cubic_41,
    fit_clpm_cubic_89,
    # RI-CLPM (within-person effect wy2 ~ wx1)
    fit_riclpm_linear,
    fit_riclpm_cubic_10,
    fit_riclpm_cubic_18,
    fit_riclpm_cubic_41,
    fit_riclpm_cubic_89,
    # DPM (y2 ~ x1)
    fit_dpm_linear,
    fit_dpm_cubic_10,
    fit_dpm_cubic_18,
    fit_dpm_cubic_41,
    fit_dpm_cubic_89,
    # CLPM on lm residuals
    fit_clpm_lm_linear,
    fit_clpm_lm_cubic_10,
    fit_clpm_lm_cubic_18,
    fit_clpm_lm_cubic_41,
    fit_clpm_lm_cubic_89,
    # CLPM on ranger residuals
    fit_clpm_rg_linear,
    fit_clpm_rg_cubic_10,
    fit_clpm_rg_cubic_18,
    fit_clpm_rg_cubic_41,
    fit_clpm_rg_cubic_89
  ),
  model_type = c(
    rep("clpm_raw",       5),
    rep("riclpm_raw",     5),
    rep("dpm_raw",        5),
    rep("clpm_lm_resid",  5),
    rep("clpm_rg_resid",  5)
  ),
  data_type = rep(c("linear", "cubic_10", "cubic_18", "cubic_41", "cubic_89"), 5),
  # which lhs and rhs correspond to "x1 -> y2"
  lhs = c(
    rep("y2",  5),  # CLPM raw
    rep("wy2", 5),  # RI-CLPM within-person: wy2 ~ wx1
    rep("y2",  5),  # DPM
    rep("y2",  5),  # CLPM on lm residuals
    rep("y2",  5)   # CLPM on ranger residuals
  ),
  rhs = c(
    rep("x1",  5),
    rep("wx1", 5),
    rep("x1",  5),
    rep("x1",  5),
    rep("x1",  5)
  )
)

# ---------------------------------------------
# extract x1 -> y2 effect from all models
# ---------------------------------------------
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

# nicely formatted tibble
effects_x1_y2_tbl <- effects_x1_y2 %>%
  arrange(model_type, data_type) %>%
  select(model_type, data_type, lhs, rhs,
         est, se, z, pvalue, ci_lower, ci_upper, std_est)

effects_x1_y2_tbl

# optional: pretty table for Quarto / R Markdown
effects_x1_y2_tbl %>%
  knitr::kable(
    digits = 3,
    caption = "Estimated cross-lagged effect x1 → y2 across models"
  )

true_value <- 0.1

bias_tbl <- effects_x1_y2_tbl %>%
  mutate(
    bias = est - true_value,
    model_type = factor(model_type),
    data_type  = factor(data_type)
  )

bias_summary <- bias_tbl %>%
  group_by(model_type, data_type) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias   = sd(bias, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(bias_tbl,
       aes(x = model_type,
           y = bias,
           fill = data_type,
           color = data_type)) +

  # boxplots: bias distribution
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.4,
    outlier.shape = NA
  ) +

  # mean ± sd per (model_type, data_type)
  geom_errorbar(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(
      x    = model_type,
      ymin = mean_bias - sd_bias,
      ymax = mean_bias + sd_bias,
      color = data_type
    ),
    width = 0.15,
    position = position_dodge(width = 0.8),
    linewidth = 0.6
  ) +

  geom_point(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(
      x = model_type,
      y = mean_bias,
      color = data_type
    ),
    size = 2.8,
    position = position_dodge(width = 0.8)
  ) +

  geom_hline(yintercept = 0, linetype = "dashed") +

  labs(
    title = "Bias of x1 → y2 across model and data types",
    x = "Model type",
    y = "Bias (estimate − 0.1)",
    color = "Data type",
    fill  = "Data type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
