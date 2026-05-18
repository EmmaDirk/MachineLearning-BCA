# ============================================================
# simple_tabpfn_clpm_once.R
# Purpose:
#   Simulate one panel dataset, residualise x/y variables using
#   TabPFN once per target, then fit a plain CLPM once.
# ============================================================

library(here)
library(reticulate)

source(here("01_scripts", "02_development", "06_thesis_v3", "00_packages.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "01_delta_sampler.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "02_delta_trajectory.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "03_simulate_panel_data.R"))
source(here("01_scripts", "02_development", "06_thesis_v3", "04_lavaan_model_string_builder.R"))

# ------------------------------------------------------------
# Python / TabPFN setup
# ------------------------------------------------------------
py_require(c("tabpfn"))

if (!py_module_available("tabpfn")) {
  stop("Python module 'tabpfn' is not available.")
}

tabpfn <- import("tabpfn")

# ------------------------------------------------------------
# settings
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# helper: confounder columns
# ------------------------------------------------------------
get_c_cols_main_only <- function(df, k = NULL, c_prefix = "c") {
  if (is.null(k)) {
    c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value = TRUE)
  } else {
    c_cols <- paste0(c_prefix, seq_len(k))
  }

  missing_c <- setdiff(c_cols, names(df))
  if (length(missing_c) > 0) {
    stop("Missing confounders: ", paste(missing_c, collapse = ", "))
  }

  c_cols
}

# ------------------------------------------------------------
# helper: TabPFN prediction for one target
# ------------------------------------------------------------
predict_tabpfn_once <- function(df, y_col, c_cols, tabpfn_module) {
  X <- as.matrix(df[, c_cols, drop = FALSE])
  y <- df[[y_col]]

  model <- tabpfn_module$TabPFNRegressor()
  model$fit(X, y)

  as.numeric(model$predict(X))
}

# ------------------------------------------------------------
# helper: residualise all x/y variables using TabPFN
# ------------------------------------------------------------
residualise_panel_tabpfn_once <- function(df,
                                          k = NULL,
                                          x_prefix = "x",
                                          y_prefix = "y",
                                          c_prefix = "c",
                                          tabpfn_module) {
  df <- as.data.frame(df)

  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value = TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value = TRUE)
  c_cols <- get_c_cols_main_only(df, k = k, c_prefix = c_prefix)

  for (v in c(x_cols, y_cols)) {
    pred <- predict_tabpfn_once(
      df = df,
      y_col = v,
      c_cols = c_cols,
      tabpfn_module = tabpfn_module
    )

    df[[v]] <- df[[v]] - pred
  }

  df
}

# ------------------------------------------------------------
# helper: fit plain CLPM
# ------------------------------------------------------------
fit_plain_clpm <- function(df, T) {
  model_string <- build_clpm(T = T, k = 0, confounder_order = 0, exclude = NULL)

  lavaan::lavaan(
    model = model_string,
    data = as.data.frame(df),
    estimator = "ML",
    warn = FALSE
  )
}

# ------------------------------------------------------------
# helper: extract lagged estimates
# ------------------------------------------------------------
extract_clpm_lags_simple <- function(fit, T) {
  out <- data.frame(
    t = 2:T,
    ARX = NA_real_,
    ARY = NA_real_,
    CXY = NA_real_,
    CYX = NA_real_
  )

  pe <- lavaan::parameterEstimates(fit)

  grab <- function(lhs, rhs) {
    ix <- which(pe$op == "~" & pe$lhs == lhs & pe$rhs == rhs)
    if (length(ix) == 0) return(NA_real_)
    pe$est[ix[1]]
  }

  for (tt in 2:T) {
    out$ARX[out$t == tt] <- grab(paste0("x", tt), paste0("x", tt - 1))
    out$ARY[out$t == tt] <- grab(paste0("y", tt), paste0("y", tt - 1))
    out$CXY[out$t == tt] <- grab(paste0("y", tt), paste0("x", tt - 1))
    out$CYX[out$t == tt] <- grab(paste0("x", tt), paste0("y", tt - 1))
  }

  out
}

# ------------------------------------------------------------
# 1. simulate one dataset
# ------------------------------------------------------------
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

cat("Simulated data:\n")
print(head(df))

# ------------------------------------------------------------
# 2. one-off TabPFN warm-up (optional but often useful)
# ------------------------------------------------------------
warm_df <- data.frame(
  c1 = rnorm(40),
  c2 = rnorm(40),
  c3 = rnorm(40)
)
warm_y <- 0.5 * warm_df$c1 - 0.3 * warm_df$c2 + rnorm(40, sd = 0.5)

warm_model <- tabpfn$TabPFNRegressor()
warm_model$fit(as.matrix(warm_df), warm_y)
invisible(warm_model$predict(as.matrix(warm_df)))

rm(warm_df, warm_y, warm_model)

# ------------------------------------------------------------
# 3. residualise x/y using TabPFN
# ------------------------------------------------------------
df_resid <- residualise_panel_tabpfn_once(
  df = df,
  k = k,
  tabpfn_module = tabpfn
)

cat("\nResidualised data:\n")
print(head(df_resid))

# ------------------------------------------------------------
# 4. fit plain CLPM once
# ------------------------------------------------------------
fit <- fit_plain_clpm(df = df_resid, T = T)

cat("\nConverged:\n")
print(lavaan::lavInspect(fit, "converged"))

cat("\nFit measures:\n")
print(lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "bic")))

cat("\nLagged estimates:\n")
print(extract_clpm_lags_simple(fit, T))
