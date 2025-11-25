### ================================================================
###  SIMULATION WITH TIME-VARYING (RANDOM) EFFECTS OF BASELINE CONFOUNDERS
###  Methods: RI-CLPM, BCA-CLPM, RF-BCA-CLPM, XGB-BCA-CLPM
### ================================================================

rm(list = ls())
set.seed(12345)

library(mvtnorm)
library(psych)
library(lavaan)
library(tidyverse)
library(ranger)
library(parallel)
library(xgboost)

### ----------------------------------------------------------------
### Multithreading settings
### ----------------------------------------------------------------

n_threads <- max(1, detectCores() - 1)   # use all cores minus one
cat("Using", n_threads, "threads for ranger and xgboost.\n")

### ----------------------------------------------------------------
### Helper: sample gamma_x and gamma_y
### ----------------------------------------------------------------

sample_gamma <- function(
    k_lin,
    k_nonlin,
    s0,
    eta0,
    min_abs = 0.01,
    max_abs = 0.30,
    max_tries = 1e6
) {

  k <- k_lin + k_nonlin

  LX <- (1 - eta0) * s0 / 2
  LY <- (1 - eta0) * s0 / 2
  NX <- eta0 * s0 / 2
  NY <- eta0 * s0 / 2

  for (try in 1:max_tries) {

    u_x_lin  <- rnorm(k_lin)
    u_y_lin  <- rnorm(k_lin)
    u_x_nlin <- rnorm(k_nonlin)
    u_y_nlin <- rnorm(k_nonlin)

    u_x_lin  <- u_x_lin  / sqrt(sum(u_x_lin^2))
    u_y_lin  <- u_y_lin  / sqrt(sum(u_y_lin^2))
    u_x_nlin <- u_x_nlin / sqrt(sum(u_x_nlin^2))
    u_y_nlin <- u_y_nlin / sqrt(sum(u_y_nlin^2))

    gamma_x_lin  <- sqrt(LX) * u_x_lin
    gamma_y_lin  <- sqrt(LY) * u_y_lin
    gamma_x_nlin <- sqrt(NX) * u_x_nlin
    gamma_y_nlin <- sqrt(NY) * u_y_nlin

    gamma_x <- c(gamma_x_lin,  gamma_x_nlin)
    gamma_y <- c(gamma_y_lin,  gamma_y_nlin)

    if (
      all(abs(gamma_x) >= min_abs & abs(gamma_x) <= max_abs &
          abs(gamma_y) >= min_abs & abs(gamma_y) <= max_abs)
    ) {
      return(list(
        gamma_x = gamma_x,
        gamma_y = gamma_y
      ))
    }
  }

  stop("Failed to find a valid sample.")
}

### ----------------------------------------------------------------
### Simulation constants
### ----------------------------------------------------------------

N <- 100000      # reduced for speed
T <- 5

ax <- 0.25
ay <- 0.25
bx <- 0.10
by <- 0.10
rho <- 0.30

A <- matrix(c(ax, by, 
              bx, ay), 
            nrow = 2, byrow = TRUE)

### ----------------------------------------------------------------
### Generate baseline confounders: linear and nonlinear
### ----------------------------------------------------------------

k_lin <- 3
k_nonlin <- 3
k <- k_lin + k_nonlin

set.seed(42)
sampled <- sample_gamma(k_lin, k_nonlin, s0 = 0.5, eta0 = 0.8)

gamma_x <- sampled$gamma_x
gamma_y <- sampled$gamma_y

gamma_x_lin  <- gamma_x[1:3]
gamma_x_tanh <- gamma_x[4:6]
gamma_y_lin  <- gamma_y[1:3]
gamma_y_tanh <- gamma_y[4:6]

# base linear confounders
U_lin <- rmvnorm(N, mean = rep(0, k_lin), sigma = diag(k_lin))
colnames(U_lin) <- paste0("c", 1:k_lin)

# tanh non-linear confounders
a_tanh <- 2
U_tanh <- apply(U_lin, 2, function(col) {
  z_raw <- tanh(a_tanh * col)
  z_perp <- resid(lm(z_raw ~ col))
  scale(z_perp, center = TRUE, scale = TRUE)[,1]
})

colnames(U_tanh) <- paste0("c", 1:k_lin, "_tanh")

U <- cbind(U_lin, U_tanh)
Psi <- cov(U)

### ----------------------------------------------------------------
### Dynamic part: compute innovation variance
### ----------------------------------------------------------------

B_avg <- rbind(gamma_x, gamma_y)
S_U <- B_avg %*% Psi %*% t(B_avg)

S_target <- diag(2)

find_c <- function(A, S_U, rho) {
  f <- function(c) {
    S_target_c <- matrix(c(1, c, c, 1), 2, 2)
    S_dyn <- S_target_c - S_U
    Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
    corr_e <- Sigma_e[1,2] / sqrt(Sigma_e[1,1] * Sigma_e[2,2])
    corr_e - rho
  }
  uniroot(f, c(-0.99, 0.99))$root
}

c_stat <- find_c(A, S_U, rho)
S_target <- matrix(c(1, c_stat, c_stat, 1), 2, 2)
S_dyn <- S_target - S_U

Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2
Sigma_e1 <- S_dyn

### ----------------------------------------------------------------
### Random wave-specific confounder effects B_t
### ----------------------------------------------------------------

B_list <- lapply(1:T, function(t) {
  m_x <- runif(1, 0.7, 1.3)
  m_y <- runif(1, 0.7, 1.3)
  rbind(
    gamma_x * m_x,
    gamma_y * m_y
  )
})

### ----------------------------------------------------------------
### Simulate data
### ----------------------------------------------------------------

df <- matrix(NA, nrow=N, ncol=2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  colnames(U))
df[, (2*T+1):(2*T+k)] <- U

# t = 1
Ddyn <- rmvnorm(N, mean=c(0,0), sigma=Sigma_e1)
B1 <- B_list[[1]]
obs1 <- Ddyn + U %*% t(B1)
df[,1]    <- obs1[,1]
df[,1+T]  <- obs1[,2]

# t = 2..T
for(i in 2:T){
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma=Sigma_e)
  Bi <- B_list[[i]]
  obs <- Ddyn + U %*% t(Bi)
  df[,i]   <- obs[,1]
  df[,i+T] <- obs[,2]
}

df <- as.data.frame(df)

### ----------------------------------------------------------------
### Fit the RI-CLPM
### ----------------------------------------------------------------

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

fit_riclpm <- lavaan::sem(
  model = model_riclpm,
  data  = df,
  meanstructure = TRUE,
  estimator     = "ML",
  missing       = "FIML",
  fixed.x       = FALSE
)

summary(fit_riclpm, standardized=TRUE, fit.measures=TRUE)

### ----------------------------------------------------------
### CORRECTED BCA RESIDUALIZATION (X and Y at every wave)
### ----------------------------------------------------------

for(t in 1:T){
  xt <- paste0("x",t)
  yt <- paste0("y",t)

  df[[paste0(xt,"_res")]] <- lm(df[[xt]] ~ df$c1 + df$c2 + df$c3)$residuals
  df[[paste0(yt,"_res")]] <- lm(df[[yt]] ~ df$c1 + df$c2 + df$c3)$residuals
}

df_res <- df %>%
  transmute(
    x1 = x1_res, x2 = x2_res, x3 = x3_res, x4 = x4_res, x5 = x5_res,
    y1 = y1_res, y2 = y2_res, y3 = y3_res, y4 = y4_res, y5 = y5_res
  )

### ----------------------------------------------------------
### CLPM ON RESIDUALS (BCA-SEM)
### ----------------------------------------------------------

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

  x1 ~ 1; x2 ~ 1; x3 ~ 1; x4 ~ 1; x5 ~ 1
  y1 ~ 1; y2 ~ 1; y3 ~ 1; y4 ~ 1; y5 ~ 1
"

fit_clpm <- lavaan::sem(
  model = model_clpm,
  data  = df_res,
  meanstructure = TRUE,
  estimator     = "ML"
)

summary(fit_clpm, standardized = TRUE, fit.measures = TRUE)

### ================================================================
### RF-BCA: BCA SEM with random forests — MULTITHREADED & FASTER
### ================================================================

## Helper to fit an RF with common settings
fit_rf_fast <- function(formula, data) {
  ranger(
    formula         = formula,
    data            = data,
    num.trees       = 1000,
    mtry            = 2,
    max.depth       = 8,
    splitrule       = "extratrees",
    honesty         = FALSE,
    sample.fraction = 0.7,
    importance      = "none",
    num.threads     = n_threads,
    verbose         = FALSE
  )
}

## Outcomes for RF
outcomes <- c(paste0("x", 1:5), paste0("y", 1:5))

cat("\nFitting RF residualizers...\n")
pb_rf <- txtProgressBar(min = 0, max = length(outcomes), style = 3)
step_rf <- 0

for (var in outcomes) {
  fml   <- as.formula(paste0(var, " ~ c1 + c2 + c3"))
  rf_mod <- fit_rf_fast(fml, df)
  preds  <- predict(rf_mod, data = df)$predictions
  df[[paste0(var, "_res_rf")]] <- df[[var]] - preds

  step_rf <- step_rf + 1
  setTxtProgressBar(pb_rf, step_rf)
}
close(pb_rf)
cat("\nRF residualization done.\n")

## Build dataset of RF residuals only for CLPM

df_res_rf <- df %>%
  transmute(
    x1 = x1_res_rf,
    x2 = x2_res_rf,
    x3 = x3_res_rf,
    x4 = x4_res_rf,
    x5 = x5_res_rf,
    y1 = y1_res_rf,
    y2 = y2_res_rf,
    y3 = y3_res_rf,
    y4 = y4_res_rf,
    y5 = y5_res_rf
  )

## CLPM on RF residuals

fit_clpm_rf <- lavaan::sem(
  model         = model_clpm,
  data          = df_res_rf,
  meanstructure = TRUE,
  estimator     = "ML"
)

cat("\n\n================ RF-BCA CLPM RESULTS ================\n")
summary(fit_clpm_rf, standardized = TRUE, fit.measures = TRUE)

###########################################################################
### XGB-BCA: XGBOOST RESIDUALIZATION (FAST)
###########################################################################

cat("\nFitting XGBoost residualizers...\n")

# Convert confounders once (matrix form)
Xmat <- as.matrix(df[, c("c1", "c2", "c3")])

# Fast XGBoost settings
xgb_params <- list(
  booster = "gbtree",
  eta = 0.10,
  max_depth = 4,
  subsample = 0.7,
  colsample_bytree = 1.0,
  objective = "reg:squarederror",
  nthread = n_threads
)

nrounds <- 80   # fast + effective

fit_xgb <- function(y) {
  xgb.train(
    params  = xgb_params,
    data    = xgb.DMatrix(Xmat, label = y),
    nrounds = nrounds,
    verbose = 0
  )
}

mods_xgb <- list()
pb_xgb <- txtProgressBar(min = 0, max = length(outcomes), style = 3)
step_xgb <- 0

for (var in outcomes) {
  mods_xgb[[var]] <- fit_xgb(df[[var]])
  df[[paste0(var, "_res_xgb")]] <- df[[var]] -
    predict(mods_xgb[[var]], Xmat)

  step_xgb <- step_xgb + 1
  setTxtProgressBar(pb_xgb, step_xgb)
}
close(pb_xgb)
cat("\nXGBoost residualization done.\n")

# Build pure residual dataset for XGB-BCA
df_res_xgb <- df %>%
  transmute(
    x1 = x1_res_xgb, x2 = x2_res_xgb, x3 = x3_res_xgb,
    x4 = x4_res_xgb, x5 = x5_res_xgb,
    y1 = y1_res_xgb, y2 = y2_res_xgb, y3 = y3_res_xgb,
    y4 = y4_res_xgb, y5 = y5_res_xgb
  )

###########################################################################
### CLPM ON XGB RESIDUALS
###########################################################################

fit_clpm_xgb <- lavaan::sem(
  model = model_clpm,
  data  = df_res_xgb,
  meanstructure = TRUE,
  estimator = "ML"
)

cat("\n\n================ XGB-BCA CLPM RESULTS ================\n")
summary(fit_clpm_xgb, standardized = TRUE, fit.measures = TRUE)

### ============================================================
### EXTRACT EFFECTS FOR ALL FOUR METHODS
### ============================================================

extract_effects <- function(fit, prefix_wx = "wx", prefix_wy = "wy",
                            prefix_x = "x", prefix_y = "y",
                            standardized = TRUE,
                            method_name) {

  pe <- parameterEstimates(fit, standardized = standardized)

  # decide whether model uses wx/wy (RI-CLPM) or plain x/y (CLPM)
  if (any(grepl("^wx", pe$lhs))) {
    lhs_x <- paste0(prefix_wx, 2:5)
    lhs_y <- paste0(prefix_wy, 2:5)
    rhs_x <- paste0(prefix_wx, 1:4)
    rhs_y <- paste0(prefix_wy, 1:4)
  } else {
    lhs_x <- paste0(prefix_x, 2:5)
    lhs_y <- paste0(prefix_y, 2:5)
    rhs_x <- paste0(prefix_x, 1:4)
    rhs_y <- paste0(prefix_y, 1:4)
  }

  out <- pe %>%
    dplyr::filter(
      op == "~",
      lhs %in% c(lhs_x, lhs_y),
      rhs %in% c(rhs_x, rhs_y)
    ) %>%
    mutate(
      effect = dplyr::case_when(
        lhs %in% lhs_x & rhs %in% rhs_x ~ "x_lag_on_x",
        lhs %in% lhs_y & rhs %in% rhs_y ~ "y_lag_on_y",
        lhs %in% lhs_x & rhs %in% rhs_y ~ "y_lag_on_x",
        lhs %in% lhs_y & rhs %in% rhs_x ~ "x_lag_on_y"
      )
    ) %>%
    group_by(effect) %>%
    summarise(
      estimate = mean(std.all),
      se       = mean(se),
      ci_low   = estimate - 1.96 * se,
      ci_high  = estimate + 1.96 * se,
      method   = method_name,
      .groups  = "drop"
    )

  return(out)
}

# extract for each model
tab_ri   <- extract_effects(fit_riclpm,   method_name = "RI-CLPM")
tab_bca  <- extract_effects(fit_clpm,     method_name = "BCA-CLPM")
tab_rf   <- extract_effects(fit_clpm_rf,  method_name = "RF-BCA-CLPM")
tab_xgb  <- extract_effects(fit_clpm_xgb, method_name = "XGB-BCA-CLPM")

# truth
truth_tbl <- tibble(
  effect = c("x_lag_on_x", "x_lag_on_y", "y_lag_on_x", "y_lag_on_y"),
  true   = c(ax,           bx,           by,           ay)
)

# combine everything
effects_all <- bind_rows(tab_ri, tab_bca, tab_rf, tab_xgb) %>%
  left_join(truth_tbl, by = "effect")

print(effects_all)

### ============================================================
### BOXPLOT WITH TRUTH LINE FOR ALL FOUR METHODS
### ============================================================

# duplicate rows 10× to make boxplots visible
effects_plot <- effects_all[rep(1:nrow(effects_all), each = 10), ]

ggplot(effects_plot, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA) +
  geom_point(
    data = effects_all,
    aes(x = method, y = estimate),
    size = 2
  ) +
  geom_errorbar(
    data = effects_all,
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_hline(
    data = effects_all,
    aes(yintercept = true),
    linetype = "dashed",
    color = "black"
  ) +
  facet_wrap(~effect, scales = "free_y") +
  labs(
    x = "",
    y = "Estimated effect",
    title = "Cross-lagged effects: RI-CLPM vs BCA-CLPM vs RF-BCA-CLPM vs XGB-BCA-CLPM"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
