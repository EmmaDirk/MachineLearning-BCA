# ============================================================
# simple_elnet_clpm_once_minimal.R
# Purpose:
#   Simulate one panel dataset, residualise x/y variables using
#   elastic net once per target, then fit a plain CLPM once.
# ============================================================

library(here)
library(glmnet)

source(here("01_scripts", "02_development", "06_thesis_v3", "00_packages.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "01_delta_sampler.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "02_delta_trajectory.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "04_lavaan_model_string_builder.R"))

seed <- 20260402L

N <- 300L
T <- 3L
k <- 3L

Phi <- matrix(c(
  0.3, 0.1,
  0.1, 0.3
), 2, 2, byrow = TRUE)

Sigma <- matrix(c(
  1.0, 0.3,
  0.3, 1.0
), 2, 2, byrow = TRUE)

Omega11 <- diag(k)

R2_total <- 0.20
R2_interaction <- 0.00
include_2way <- FALSE
include_3way <- FALSE
burn_in <- 0L

alpha_elnet <- 0.50

d1_obj <- sample_delta_1(
  k = k,
  Omega11 = Omega11,
  R2_total = R2_total,
  R2_interaction = R2_interaction,
  include_2way = include_2way,
  include_3way = include_3way
)

Delta1 <- d1_obj$Delta

Delta_list <- generate_Delta_constant(
  Delta1 = Delta1,
  T = T,
  burn_in = burn_in
)

df <- simulate_panel_data(
  N = N,
  T = T,
  Phi = Phi,
  Delta_list = Delta_list,
  Omega11 = Omega11,
  Sigma = Sigma,
  burn_in = burn_in,
  seed = seed
)

x_cols <- grep("^x\\d+$", names(df), value = TRUE)
y_cols <- grep("^y\\d+$", names(df), value = TRUE)
c_cols <- paste0("c", seq_len(k))

for (v in c(x_cols, y_cols)) {
  X <- as.matrix(df[, c_cols, drop = FALSE])
  y <- df[[v]]

  set.seed(seed)

  cv_fit <- glmnet::cv.glmnet(
    x = X,
    y = y,
    alpha = alpha_elnet,
    family = "gaussian",
    standardize = TRUE,
    intercept = TRUE,
    nfolds = 5
  )

  pred <- as.numeric(predict(cv_fit, newx = X, s = "lambda.min"))
  df[[v]] <- y - pred
}

model_string <- build_clpm(
  T = T,
  k = 0,
  confounder_order = 0,
  exclude = NULL
)

fit <- lavaan::lavaan(
  model = model_string,
  data = as.data.frame(df),
  estimator = "ML",
  warn = FALSE
)

pe <- lavaan::parameterEstimates(fit)

lagged_estimates <- data.frame(
  t = 2:T,
  ARX = NA_real_,
  ARY = NA_real_,
  CXY = NA_real_,
  CYX = NA_real_
)

for (tt in 2:T) {
  ix_arx <- which(pe$op == "~" & pe$lhs == paste0("x", tt) & pe$rhs == paste0("x", tt - 1))
  ix_ary <- which(pe$op == "~" & pe$lhs == paste0("y", tt) & pe$rhs == paste0("y", tt - 1))
  ix_cxy <- which(pe$op == "~" & pe$lhs == paste0("y", tt) & pe$rhs == paste0("x", tt - 1))
  ix_cyx <- which(pe$op == "~" & pe$lhs == paste0("x", tt) & pe$rhs == paste0("y", tt - 1))

  if (length(ix_arx) > 0) lagged_estimates$ARX[lagged_estimates$t == tt] <- pe$est[ix_arx[1]]
  if (length(ix_ary) > 0) lagged_estimates$ARY[lagged_estimates$t == tt] <- pe$est[ix_ary[1]]
  if (length(ix_cxy) > 0) lagged_estimates$CXY[lagged_estimates$t == tt] <- pe$est[ix_cxy[1]]
  if (length(ix_cyx) > 0) lagged_estimates$CYX[lagged_estimates$t == tt] <- pe$est[ix_cyx[1]]
}

head(df)
lavaan::lavInspect(fit, "converged")
lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "bic"))
lagged_estimates
