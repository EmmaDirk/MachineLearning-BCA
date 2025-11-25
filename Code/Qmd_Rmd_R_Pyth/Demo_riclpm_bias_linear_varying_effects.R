# ================================================================
# Linear panel-data simulation with time-varying confounder effects
# RI-CLPM estimation  +  Corrected BCA-CLPM
# ================================================================

rm(list = ls())
set.seed(427)

library(mvtnorm)
library(psych)
library(lavaan)
library(tidyverse)

### ----------------------------------------------------------
### MODEL PARAMETERS
### ----------------------------------------------------------

N <- 1000000
T <- 5

ax <- 0.25
ay <- 0.25
bx <- 0.10
by <- 0.10
rho <- 0.30

A <- matrix(c(ax, by,
              bx, ay), 2, 2, byrow = TRUE)

### ----------------------------------------------------------
### LINEAR CONFOUNDERS
### ----------------------------------------------------------

k <- 3
tau2 <- rep(1, k)
Psi <- diag(tau2)

gamma_x_base <- c(0.40, 0.30, 0.20)
gamma_y_base <- c(0.40, 0.30, 0.20)

# ---------------------------------------------------------------
# Time-varying confounder effects
# ---------------------------------------------------------------
set.seed(777)

gamma_x_list <- lapply(1:T, function(t) gamma_x_base * runif(1, 0.7, 1.3))
gamma_y_list <- lapply(1:T, function(t) gamma_y_base * runif(1, 0.7, 1.3))

B_list <- lapply(1:T, function(t) {
  rbind(gamma_x_list[[t]], gamma_y_list[[t]])
})

### ----------------------------------------------------------
### STATIONARITY CALCULATIONS (using mean B)
### ----------------------------------------------------------

B_avg <- rbind(gamma_x_base, gamma_y_base)
S_U <- B_avg %*% Psi %*% t(B_avg)

S_target <- diag(2)

find_c <- function(A, S_U, rho) {
  f <- function(c) {
    S_target_c <- matrix(c(1, c,
                           c, 1), 2, 2)
    S_dyn_c <- S_target_c - S_U
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    corr_e <- Sigma_e_c[1,2] / sqrt(Sigma_e_c[1,1] *
                                    Sigma_e_c[2,2])
    corr_e - rho
  }
  uniroot(f, interval = c(-0.99, 0.99))$root
}

c_stat <- find_c(A, S_U, rho)

S_target <- matrix(c(1, c_stat,
                     c_stat, 1), 2, 2)

S_dyn <- S_target - S_U
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2

Sigma_e1 <- S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

### ----------------------------------------------------------
### SIMULATE CONFOUNDERS
### ----------------------------------------------------------

U <- rmvnorm(N, rep(0, k), Psi)

### ----------------------------------------------------------
### SIMULATE PANEL DATA X1..X5, Y1..Y5
### ----------------------------------------------------------

df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x",1:T),
                  paste0("y",1:T),
                  paste0("c",1:k))
df[, (2*T+1):(2*T+k)] <- U

# t = 1
Ddyn <- rmvnorm(N, c(0,0), Sigma_e1)
B1 <- B_list[[1]]
obs1 <- Ddyn + U %*% t(B1)

df[,1]     <- obs1[,1]
df[,1+T]   <- obs1[,2]

# t = 2..T
for(i in 2:T){
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e)
  Bi <- B_list[[i]]
  obs <- Ddyn + U %*% t(Bi)

  df[, i]     <- obs[,1]
  df[, i + T] <- obs[,2]
}

df <- as.data.frame(df)

### ----------------------------------------------------------
### Basic checks
### ----------------------------------------------------------

print(round(apply(df, 2, var), 3)[1:10])
print(round(cov(df)[1:10, 1:10], 3))

### ----------------------------------------------------------
### RI-CLPM
### ----------------------------------------------------------

model_riclpm <- "
  # Random intercepts
  rix =~ x1 + x2 + x3 + x4 + x5
  riy =~ y1 + y2 + y3 + y4 + y5
  rix ~~ rix
  riy ~~ riy
  rix ~~ riy

  # Observed residual variances forced to RI + within
  x1 ~~ 0*x1; x2 ~~ 0*x2; x3 ~~ 0*x3; x4 ~~ 0*x4; x5 ~~ 0*x5
  y1 ~~ 0*y1; y2 ~~ 0*y2; y3 ~~ 0*y3; y4 ~~ 0*y4; y5 ~~ 0*y5

  # Within-person factors
  wx1 =~ 1*x1; wx2 =~ 1*x2; wx3 =~ 1*x3; wx4 =~ 1*x4; wx5 =~ 1*x5
  wy1 =~ 1*y1; wy2 =~ 1*y2; wy3 =~ 1*y3; wy4 =~ 1*y4; wy5 =~ 1*y5

  # RI orthogonality
  rix ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
  rix ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5
  riy ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
  riy ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5

  # Within residual variances
  wx1 ~~ wx1; wx2 ~~ wx2; wx3 ~~ wx3; wx4 ~~ wx4; wx5 ~~ wx5
  wy1 ~~ wy1; wy2 ~~ wy2; wy3 ~~ wy3; wy4 ~~ wy4; wy5 ~~ wy5
  wy1 ~~ wx1; wy2 ~~ wx2; wy3 ~~ wx3; wy4 ~~ wx4; wy5 ~~ wx5

  # Cross-lagged structure
  wx2 ~ wx1 + wy1
  wy2 ~ wx1 + wy1
  wx3 ~ wx2 + wy2
  wy3 ~ wx2 + wy2
  wx4 ~ wx3 + wy3
  wy4 ~ wx3 + wy3
  wx5 ~ wx4 + wy4
  wy5 ~ wx4 + wy4

  # Means
  x1 + x2 + x3 + x4 + x5 ~ 1
  y1 + y2 + y3 + y4 + y5 ~ 1
"

fit_riclpm <- lavaan::sem(
  model = model_riclpm,
  data = df,
  meanstructure = TRUE,
  estimator = "ML"
)

### ----------------------------------------------------------
### CORRECTED BCA RESIDUALIZATION (X and Y at every wave)
### ----------------------------------------------------------

# residualize X_t ~ C, Y_t ~ C
for(t in 1:T){
  xt <- paste0("x",t)
  yt <- paste0("y",t)

  df[[paste0(xt,"_res")]] <- lm(df[[xt]] ~ df$c1 + df$c2 + df$c3)$residuals
  df[[paste0(yt,"_res")]] <- lm(df[[yt]] ~ df$c1 + df$c2 + df$c3)$residuals
}

# build pure-residual dataset
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
  data = df_res,
  meanstructure = TRUE,
  estimator = "ML"
)

summary(fit_riclpm, standardized = TRUE, fit.measures = TRUE)

summary(fit_clpm, standardized = TRUE, fit.measures = TRUE)

# ------------------------------
# TRUE EFFECTS
# ------------------------------
truth <- tibble(
  effect = c("x_lag_on_x", "y_lag_on_y", "y_lag_on_x", "x_lag_on_y"),
  true   = c(ax, ay, by, bx)
)

# ------------------------------
# RI-CLPM extraction (std)
# ------------------------------
pe_ric <- parameterEstimates(fit_riclpm, standardized = TRUE) %>%
  filter(op == "~",
         lhs %in% paste0("wx", 2:5) | lhs %in% paste0("wy", 2:5),
         rhs %in% paste0("wx", 1:4) | rhs %in% paste0("wy", 1:4)) %>%
  mutate(
    effect = case_when(
      lhs %in% paste0("wx", 2:5) & rhs %in% paste0("wx", 1:4) ~ "x_lag_on_x",
      lhs %in% paste0("wy", 2:5) & rhs %in% paste0("wy", 1:4) ~ "y_lag_on_y",
      lhs %in% paste0("wx", 2:5) & rhs %in% paste0("wy", 1:4) ~ "y_lag_on_x",
      lhs %in% paste0("wy", 2:5) & rhs %in% paste0("wx", 1:4) ~ "x_lag_on_y"
    )
  ) %>%
  group_by(effect) %>%
  summarise(
    estimate = mean(std.all),
    se       = mean(se),
    ci_low   = estimate - 1.96 * se,
    ci_high  = estimate + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(method = "RI-CLPM")

# ------------------------------
# BCA-CLPM extraction (std)
# ------------------------------
pe_clp <- parameterEstimates(fit_clpm, standardized = TRUE) %>%
  filter(op == "~",
         lhs %in% paste0("x", 2:5) | lhs %in% paste0("y", 2:5),
         rhs %in% paste0("x", 1:4) | rhs %in% paste0("y", 1:4)) %>%
  mutate(
    effect = case_when(
      lhs %in% paste0("x", 2:5) & rhs %in% paste0("x", 1:4) ~ "x_lag_on_x",
      lhs %in% paste0("y", 2:5) & rhs %in% paste0("y", 1:4) ~ "y_lag_on_y",
      lhs %in% paste0("x", 2:5) & rhs %in% paste0("y", 1:4) ~ "y_lag_on_x",
      lhs %in% paste0("y", 2:5) & rhs %in% paste0("x", 1:4) ~ "x_lag_on_y"
    )
  ) %>%
  group_by(effect) %>%
  summarise(
    estimate = mean(std.all),
    se       = mean(se),
    ci_low   = estimate - 1.96 * se,
    ci_high  = estimate + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(method = "BCA-CLPM")

# ------------------------------
# Combined effects table
# ------------------------------
effects_table <- bind_rows(pe_ric, pe_clp) %>%
  left_join(truth, by = "effect")

effects_table

# Aesthetics: boxplot requires multiple rows per method, but we have one row each.
# We intentionally "fake" a narrow box by duplicating each row.
plot_data <- effects_table %>%
  group_by(effect, method) %>%
  summarize(
    estimate = estimate,
    ci_low = ci_low,
    ci_high = ci_high,
    true = true
  ) %>%
  ungroup()

# Fake 10 copies to create a visible box (same value each time)
plot_data_long <- plot_data[rep(1:nrow(plot_data), each = 10), ]

ggplot(plot_data_long, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA) +
  geom_point(data = plot_data, aes(y = estimate), size = 2) +
  geom_errorbar(
    data = plot_data,
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_hline(
    data = plot_data,
    aes(yintercept = true),
    linetype = "dashed",
    color = "black"
  ) +
  facet_wrap(~effect, scales = "free_y") +
  labs(
    x = "",
    y = "Estimated effect",
    title = "RI-CLPM vs BCA-CLPM — Cross-lagged Effects"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
