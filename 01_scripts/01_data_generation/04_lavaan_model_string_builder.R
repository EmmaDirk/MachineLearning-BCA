# =================================================================================================
#
# These functions build the lavaan model strings for the different models we want to fit.
# This needs to be done to adapt the models to the number of time points T.
# As such, we use text manipulation to create the model strings dynamically.
#
# The functions build the following models:
# - 1) CLPM where we can optionally control for confounders directly
# - 2) RI-CLPM, optionally freeing loadings
# - 3) DPM, optionally freeing loadings
# =================================================================================================

# ---- observed confounder terms ------------------------------------------------------------------

# confounder_order controls how observed confounders enter the model:
# 0 = no confounders
# 1 = main effects only
# 2 = main effects + all 2-way interactions among the confounders
# 3 = main effects + all 2-way + all 3-way interactions among the confounders
#
# exclude optionally removes specific confounder terms, such as "C1" or "Z1.2".

build_observed_confounder_terms <- function(k = 0, confounder_order = 0, exclude = NULL) {

  # Build analyst-side confounder terms.
  # The names follow the Appendix A blocks C, Z, and T after lavaan-safe dot conversion.

  # The C block contains the observed base confounders C1, C2, ..., Ck.

  C_names <- if (k > 0) paste0("C", 1:k) else character(0)

  # Excluded terms are removed after they have been converted to the same lavaan-safe form.

  remove_excluded <- function(x) {
    if (is.null(exclude)) return(x)
    setdiff(x, exclude)
  }

  # conf_terms is built in the same conceptual order as M = (C', Z', T')'.

  conf_terms <- character(0)

  # Order 1 includes only the base-confounder block C.

  if (confounder_order >= 1 && k > 0) {
    main_terms <- remove_excluded(C_names)
    conf_terms <- c(conf_terms, main_terms)
  }

  # Order 2 additionally includes standardized two-way interaction terms Z.

  if (confounder_order >= 2 && k >= 2) {
    two_way <- combn(seq_len(k), 2, FUN = function(ix) paste0("Z", paste(ix, collapse = ".")))
    two_way <- remove_excluded(two_way)
    conf_terms <- c(conf_terms, two_way)
  }

  # Order 3 additionally includes standardized three-way interaction terms T.

  if (confounder_order >= 3 && k >= 3) {
    three_way <- combn(seq_len(k), 3, FUN = function(ix) paste0("T", paste(ix, collapse = ".")))
    three_way <- remove_excluded(three_way)
    conf_terms <- c(conf_terms, three_way)
  }

  conf_terms
}


# ---- CLPM ----------------------------------------------------------------------------------------

# If confounders are included, they enter:
# - wave 1 directly: x1 ~ C, y1 ~ C
# - waves 2..T directly as well: x_t ~ x_{t-1} + y_{t-1} + C,
#   y_t ~ x_{t-1} + y_{t-1} + C.

build_clpm <- function(T, k = 0, confounder_order = 0, exclude = NULL) {

  # Build the observed confounder terms once.

  conf_terms <- build_observed_confounder_terms(
    k = k,
    confounder_order = confounder_order,
    exclude = exclude
  )

  # Baseline equations.
  # If confounders are included, wave 1 is also regressed on them.

  baseline_block <- character(0)

  if (length(conf_terms) > 0) {
    conf_rhs <- paste(conf_terms, collapse = " + ")
    baseline_block <- c(
      sprintf("x1 ~ %s", conf_rhs),
      sprintf("y1 ~ %s", conf_rhs)
    )
  }

  # Lagged regressions for waves 2..T.

  regress_block <- paste(
    c(
      baseline_block,
      unlist(lapply(2:T, function(t) {

        # Define lagged predictors.

        lag_x <- sprintf("x%d", t - 1)
        lag_y <- sprintf("y%d", t - 1)

        # Start with baseline predictors.

        rhs_terms <- c(lag_x, lag_y)

        # Add the same observed confounder terms at every later wave.

        if (length(conf_terms) > 0) {
          rhs_terms <- c(rhs_terms, conf_terms)
        }

        # Collapse RHS.

        rhs <- paste(rhs_terms, collapse = " + ")

        c(
          # X_t regressed on predictors.

          sprintf("x%d ~ %s", t, rhs),

          # Y_t regressed on predictors.

          sprintf("y%d ~ %s", t, rhs)
        )
      }))
    ),
    collapse = "\n"
  )

  # Residual covariances: X_t ~~ Y_t.

  resid_cov <- paste(sprintf("x%d ~~ y%d", 1:T, 1:T), collapse = "\n")

  # Residual variances.

  resid_vars <- paste(
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse = "\n"),
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse = "\n"),
    sep = "\n"
  )

  # Means.

  means_block <- paste(
    paste(paste0("x", 1:T), collapse = " + "), "~ 1\n",
    paste(paste0("y", 1:T), collapse = " + "), "~ 1\n"
  )

  # Combine all blocks.

  paste(regress_block, resid_cov, resid_vars, means_block, sep = "\n")
}


# ---- RI-CLPM -------------------------------------------------------------------------------------

# If confounders are included here, we control for them on the observed level first.
# This means we regress x1..xT and y1..yT directly on the observed confounder columns,
# before the RI / within-person decomposition is defined through the latent structure.
#
# For each observed confounder term, the corresponding coefficient is constrained
# to be equal over time within the X block and within the Y block.

build_riclpm <- function(T, free_loadings = FALSE, k = 0, confounder_order = 0, exclude = NULL) {

  if (T < 2) stop("T must be at least 2.")

  # Build the observed confounder terms once.

  conf_terms <- build_observed_confounder_terms(
    k = k,
    confounder_order = confounder_order,
    exclude = exclude
  )

  # Random intercepts.

  if (free_loadings) {
    rix_terms <- c("NA*x1", if (T >= 2) paste0("x", 2:T))
    riy_terms <- c("NA*y1", if (T >= 2) paste0("y", 2:T))

    ri_block <- paste0(
      "rix =~ ", paste(rix_terms, collapse = " + "), "\n",
      "riy =~ ", paste(riy_terms, collapse = " + "), "\n",
      "rix ~~ 1*rix\n",
      "riy ~~ 1*riy\n",
      "rix ~~ riy\n"
    )
  } else {
    rix_terms <- sprintf("1*x%d", 1:T)
    riy_terms <- sprintf("1*y%d", 1:T)

    ri_block <- paste0(
      "rix =~ ", paste(rix_terms, collapse = " + "), "\n",
      "riy =~ ", paste(riy_terms, collapse = " + "), "\n",
      "rix ~~ rix\n",
      "riy ~~ riy\n",
      "rix ~~ riy\n"
    )
  }

  # Direct observed confounder control before the latent decomposition.

  conf_block <- character(0)

  if (length(conf_terms) > 0) {

    # Constrain each observed confounder coefficient to be equal over time within X and within Y.

    conf_labels <- gsub("\\.", "_", conf_terms)
    x_rhs <- paste(sprintf("sx_%s*%s", conf_labels, conf_terms), collapse = " + ")
    y_rhs <- paste(sprintf("sy_%s*%s", conf_labels, conf_terms), collapse = " + ")

    conf_block <- c(
      sprintf("%s ~ %s", paste(paste0("x", 1:T), collapse = " + "), x_rhs),
      sprintf("%s ~ %s", paste(paste0("y", 1:T), collapse = " + "), y_rhs)
    )
  }

  # Fix observed residual variances to zero.

  resid_fix <- paste0(
    paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "; "), "\n",
    paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "; "), "\n"
  )

  # Within-person latent variables.

  within_lat <- paste0(
    paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "; "), "\n",
    paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "; "), "\n"
  )

  # Orthogonality constraints.

  orth <- paste0(
    "rix ~~ ", paste(sprintf("0*wx%d", 1:T), collapse = " + "), "\n",
    "rix ~~ ", paste(sprintf("0*wy%d", 1:T), collapse = " + "), "\n",
    "riy ~~ ", paste(sprintf("0*wx%d", 1:T), collapse = " + "), "\n",
    "riy ~~ ", paste(sprintf("0*wy%d", 1:T), collapse = " + "), "\n"
  )

  # Within-person variances.

  within_var <- paste0(
    paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse = "; "), "\n",
    paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse = "; "), "\n"
  )

  # Within-person covariances.

  within_cov <- paste0(
    paste(sprintf("wy%d ~~ wx%d", 1:T, 1:T), collapse = "; "), "\n"
  )

  # Autoregressive and cross-lagged paths.

  regress <- paste(
    unlist(lapply(2:T, function(t) {
      c(
        sprintf("wx%d ~ wx%d + wy%d", t, t - 1, t - 1),
        sprintf("wy%d ~ wx%d + wy%d", t, t - 1, t - 1)
      )
    })),
    collapse = "\n"
  )

  # Means.

  means <- if (free_loadings) {
    paste0(
      paste(sprintf("x%d ~ 1", 1:T), collapse = "; "), "\n",
      paste(sprintf("y%d ~ 1", 1:T), collapse = "; "), "\n"
    )
  } else {
    paste0(
      paste(paste0("x", 1:T), collapse = " + "), " ~ mx*1\n",
      paste(paste0("y", 1:T), collapse = " + "), " ~ my*1\n"
    )
  }

  paste(
    paste(conf_block, collapse = "\n"),
    ri_block,
    resid_fix,
    within_lat,
    orth,
    within_var,
    within_cov,
    regress,
    means,
    sep = "\n"
  )
}


# ---- DPM -----------------------------------------------------------------------------------------

# If confounders are included here, they enter on the observed level:
# - wave 1 directly: x1 ~ C, y1 ~ C
# - waves 2..T directly as well: x_t ~ x_{t-1} + y_{t-1} + C,
#   y_t ~ x_{t-1} + y_{t-1} + C.

build_dpm <- function(T, free_loadings = FALSE, k = 0, confounder_order = 0, exclude = NULL) {

  if (T < 2) stop("T must be at least 2.")

  # Build the observed confounder terms once.

  conf_terms <- build_observed_confounder_terms(
    k = k,
    confounder_order = confounder_order,
    exclude = exclude
  )

  # Accumulating factors.

  if (free_loadings) {
    fx_terms <- c("NA*x2", if (T >= 3) paste0("x", 3:T))
    fy_terms <- c("NA*y2", if (T >= 3) paste0("y", 3:T))

    latent_block <- paste0(
      "FX =~ ", paste(fx_terms, collapse = " + "), "\n",
      "FY =~ ", paste(fy_terms, collapse = " + "), "\n",
      "FX ~~ 1*FX\n",
      "FY ~~ 1*FY\n",
      "FX ~~ FY\n"
    )
  } else {
    fx_terms <- sprintf("1*x%d", 2:T)
    fy_terms <- sprintf("1*y%d", 2:T)

    latent_block <- paste0(
      "FX =~ ", paste(fx_terms, collapse = " + "), "\n",
      "FY =~ ", paste(fy_terms, collapse = " + "), "\n"
    )
  }

  # Residual covariances between factors and baseline variables.

  fx_cov_block <- "FX ~~ x1 + y1\n"
  fy_cov_block <- "FY ~~ x1 + y1\n"

  # Baseline equations.
  # If confounders are included, wave 1 is also regressed on them.

  baseline_block <- character(0)

  if (length(conf_terms) > 0) {
    conf_rhs <- paste(conf_terms, collapse = " + ")
    baseline_block <- c(
      sprintf("x1 ~ %s", conf_rhs),
      sprintf("y1 ~ %s", conf_rhs)
    )
  }

  # Autoregressive and cross-lagged paths.

  regress_block <- paste(
    c(
      baseline_block,
      sapply(2:T, function(t) {

        # Start with lagged predictors.

        rhs_terms <- c(sprintf("x%d", t - 1), sprintf("y%d", t - 1))

        # Add the same observed confounder terms at every later wave.

        if (length(conf_terms) > 0) {
          rhs_terms <- c(rhs_terms, conf_terms)
        }

        sprintf("x%d + y%d ~ %s", t, t, paste(rhs_terms, collapse = " + "))
      })
    ),
    collapse = "\n"
  )

  # Residual covariances between x_t and y_t.

  resid_cov_block <- paste(
    sprintf("x%d ~~ y%d", 1:T, 1:T),
    collapse = "\n"
  )

  # Latent variances and covariance for non-free version only.

  latent_cov_block <- if (!free_loadings) {
    paste(
      "FX ~~ FX",
      "FY ~~ FY",
      "FX ~~ FY",
      sep = "\n"
    )
  } else {
    NULL
  }

  # Residual variances.

  resid_var_block <- paste(
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse = "\n"),
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse = "\n"),
    sep = "\n"
  )

  # Means.

  means_block <- paste(
    c(sprintf("x%d ~ 1", 1:T), sprintf("y%d ~ 1", 1:T)),
    collapse = "\n"
  )

  # Combine all blocks.

  parts <- c(
    latent_block,
    fx_cov_block,
    fy_cov_block,
    regress_block,
    resid_cov_block,
    latent_cov_block,
    resid_var_block,
    means_block
  )

  paste(parts[!sapply(parts, is.null)], collapse = "\n\n")
}