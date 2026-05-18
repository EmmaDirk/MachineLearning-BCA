# =================================================================================================
#
# Appendix A notation:
#
#   eta_t = Delta_t M,
#   M     = (C', Z', T')'.
#
# C contains the standardized base confounders. Z contains standardized two-way interaction terms.
# T contains standardized three-way interaction terms. Omega denotes the covariance matrix of M.
#
# This script implements three trajectory types:
#
# 1. Constant effects.
#    Delta_t is kept fixed for every generated wave.
#
# 2. Stepwise effects.
#    Delta_t is kept fixed before a chosen observed wave and is then multiplied by one common
#    scale factor. This changes Var(eta_t) from R2_old to R2_new while preserving the complete
#    coefficient pattern.
#
# 3. Stepwise mixture effects.
#    The post-step Delta_t is constructed as a convex mixture of the initial coefficient pattern
#    and a randomly sampled alternative pattern. The mixed direction is then rescaled using the
#    Appendix A variance decomposition. This allows the rank order of the coefficients to change
#    while preserving the old main-effect versus interaction variance split.
#
# Caution about R2_old and R2_new:
#
# R2_old and R2_new refer only to the direct variance of eta_t at the wave where the corresponding
# Delta_t is applied. They do not necessarily equal the total baseline-confounder-related variance
# in X_t or Y_t after autoregressive and cross-lagged propagation. If a burn-in period is used and
# the retained waves are interpreted after the process has moved toward stationarity, the realised
# confounder-related variance in those retained waves can differ from R2_old and R2_new. These
# functions construct Delta_t trajectories only; they do not verify the propagated variance of the
# full longitudinal system.
#
# Caution about burn-in labels:
#
# burn_in is the number of unobserved Delta_t matrices prepended before the observed waves. The
# returned list keeps this convention in its positions and uses negative-style labels for clarity.
# For example, burn_in = 20 and n_waves = 5 yields labels t-19, ..., t0, t1, ..., t5. Thus, t0 is
# the final burn-in wave and t1 is the first observed wave. The number of returned matrices remains
# burn_in + n_waves, so this naming convention does not change the simulation logic.
# =================================================================================================


# ---- Shared input validation ---------------------------------------------------------------------
# These helpers check the arguments that are used by more than one trajectory function.

validate_Delta_trajectory_inputs <- function(n_waves, burn_in) {

  # n_waves is the number of observed waves requested for analysis.
  if (!is.numeric(n_waves) || length(n_waves) != 1 || is.na(n_waves) ||
      n_waves < 1 || n_waves != as.integer(n_waves)) {
    stop("n_waves must be a positive integer.")
  }

  # burn_in is the number of unobserved waves prepended to the observed trajectory.
  if (!is.numeric(burn_in) || length(burn_in) != 1 || is.na(burn_in) ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("burn_in must be a non-negative integer.")
  }
}


validate_R2 <- function(R2, name, lower_open = FALSE) {

  # R2-type inputs are scalar variance targets or variance shares.
  if (!is.numeric(R2) || length(R2) != 1 || is.na(R2)) {
    stop(paste0(name, " must be a single numeric value."))
  }

  # Some ratios require a strictly positive lower bound, whereas target variances may be zero.
  if (lower_open) {
    if (R2 <= 0 || R2 > 1) {
      stop(paste0(name, " must be in (0, 1]."))
    }
  } else {
    if (R2 < 0 || R2 > 1) {
      stop(paste0(name, " must be in [0, 1]."))
    }
  }
}


validate_lambda <- function(lambda, name) {

  # Lambda controls how strongly the post-step direction is pulled toward a random direction.
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) ||
      lambda < 0 || lambda > 1) {
    stop(paste0(name, " must be a single number in [0, 1]."))
  }
}


# ---- Burn-in construction and trajectory names ---------------------------------------------------
# These helpers prepend burn-in matrices and name the returned trajectory consistently.

prepend_burn_in <- function(Delta_observed_list, Delta_initial, burn_in) {

  # With no burn-in, the observed trajectory is already the full internal trajectory.
  if (burn_in == 0L) {
    return(Delta_observed_list)
  }

  # Burn-in waves use the initial coefficient matrix so that steps occur only in observed time.
  Delta_burn_in_list <- vector("list", burn_in)

  for (b in seq_len(burn_in)) {
    Delta_burn_in_list[[b]] <- Delta_initial
  }

  c(Delta_burn_in_list, Delta_observed_list)
}


name_Delta_trajectory <- function(Delta_list, burn_in) {

  # The returned list always contains burn_in + n_waves matrices.
  n_total <- length(Delta_list)
  n_observed <- n_total - burn_in

  # Burn-in matrices receive non-positive names, and observed matrices receive positive names.
  if (burn_in > 0L) {
    label_values <- c(seq.int(from = -burn_in + 1L, to = 0L), seq_len(n_observed))
  } else {
    label_values <- seq_len(n_observed)
  }

  names(Delta_list) <- paste0("t", label_values)

  # Only structural metadata is stored; no separate time-index column is added downstream.
  attr(Delta_list, "burn_in") <- burn_in
  attr(Delta_list, "n_total") <- n_total
  attr(Delta_list, "n_observed") <- n_observed

  Delta_list
}


# ---- Block identification for M = (C', Z', T')' --------------------------------------------------
# This helper finds the Appendix A blocks from the Delta_t column names.

identify_M_blocks <- function(Delta_initial, Omega) {

  # Delta_initial must be a matrix with one row per outcome and one column per feature in M.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  # Omega must be the covariance matrix of M in the same feature order as Delta_initial.
  if (!is.matrix(Omega)) {
    stop("Omega must be a matrix.")
  }

  # Delta_initial and Omega must refer to the same feature vector M.
  if (ncol(Delta_initial) != nrow(Omega) || nrow(Omega) != ncol(Omega)) {
    stop("Delta_initial and Omega are not conformable.")
  }

  # Feature names are required because C, Z, and T are identified by their prefixes.
  if (is.null(colnames(Delta_initial))) {
    stop("Delta_initial must have column names following the C/Z/T naming scheme.")
  }

  if (is.null(rownames(Omega)) || is.null(colnames(Omega))) {
    stop("Omega must have row and column names following the C/Z/T naming scheme.")
  }

  # Omega and Delta_initial must use the same names in the same order.
  if (!identical(colnames(Delta_initial), rownames(Omega)) ||
      !identical(colnames(Delta_initial), colnames(Omega))) {
    stop("The feature names of Delta_initial and Omega must match exactly and be in the same order.")
  }

  # C names identify base confounders, Z names identify two-way terms, and T names identify three-way terms.
  feature_names <- colnames(Delta_initial)
  idx_C <- grep("^C[0-9]+$", feature_names)
  idx_Z <- grep("^Z[0-9]+:[0-9]+$", feature_names)
  idx_T <- grep("^T[0-9]+:[0-9]+:[0-9]+$", feature_names)

  # At least one base-confounder column is required.
  if (length(idx_C) == 0L) {
    stop("No C block found. Column names should include C1, C2, ..., Ck.")
  }

  # Every feature must belong to exactly one Appendix A block.
  idx_all <- c(idx_C, idx_Z, idx_T)
  if (length(idx_all) != length(feature_names)) {
    stop("Every column name must follow the C/Z/T naming scheme.")
  }

  # The expected Appendix A order is C first, then Z, then T.
  if (!identical(idx_all, seq_along(feature_names))) {
    stop("Columns must be ordered as M = (C', Z', T')'.")
  }

  list(
    idx_C = idx_C,
    idx_Z = idx_Z,
    idx_T = idx_T,
    has_Z = length(idx_Z) > 0L,
    has_T = length(idx_T) > 0L,
    has_int = (length(idx_Z) + length(idx_T)) > 0L
  )
}


# ---- Appendix A variance decomposition -----------------------------------------------------------
# These helpers compute the direct variance of eta_t and its effective components.

qform <- function(x, Omega_block) {

  # This returns the scalar quadratic form x' Omega_block x.
  as.numeric(t(x) %*% Omega_block %*% x)
}


compute_eta_components <- function(delta_row, Omega, blocks) {

  # One row of Delta_t is split into C, Z, and T coefficient blocks.
  delta_C <- delta_row[blocks$idx_C]
  delta_Z <- if (blocks$has_Z) delta_row[blocks$idx_Z] else numeric(0)
  delta_T <- if (blocks$has_T) delta_row[blocks$idx_T] else numeric(0)

  # The covariance blocks are extracted from Omega using the same C/Z/T partition.
  Omega11 <- Omega[blocks$idx_C, blocks$idx_C, drop = FALSE]
  Omega22 <- if (blocks$has_Z) Omega[blocks$idx_Z, blocks$idx_Z, drop = FALSE] else matrix(0, 0, 0)
  Omega33 <- if (blocks$has_T) Omega[blocks$idx_T, blocks$idx_T, drop = FALSE] else matrix(0, 0, 0)
  Omega13 <- if (blocks$has_T) {
    Omega[blocks$idx_C, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, length(blocks$idx_C), 0)
  }

  # A is the C-block variance before the C/T cross term is added.
  A <- qform(delta_C, Omega11)

  # B is the Z-block variance.
  B <- if (blocks$has_Z) qform(delta_Z, Omega22) else 0

  # G is the T-block variance.
  G <- if (blocks$has_T) qform(delta_T, Omega33) else 0

  # D is the C/T covariance term.
  D <- if (blocks$has_T) as.numeric(t(delta_C) %*% Omega13 %*% delta_T) else 0

  # Appendix A splits the total C/T cross term symmetrically between both components.
  Vlin <- A + D
  Vint <- B + G + D

  # The total direct variance of eta_t is the sum of the effective components.
  Vtotal <- Vlin + Vint

  list(A = A, B = B, G = G, D = D, Vlin = Vlin, Vint = Vint, Vtotal = Vtotal)
}


# ---- Appendix A scale factors --------------------------------------------------------------------
# This helper solves for sL and s after a new coefficient direction has been chosen.

solve_appendix_scales <- function(A, B, G, D, Vlin_star, Vint_star, tol = 1e-10) {

  # sL scales the main-effect direction b_L, and s scales the interaction directions b_2 and b_3.
  # The target equations are Vlin_star = sL^2 A + sL s D and Vint_star = s^2(B + G) + sL s D.

  # If the target total is zero, both coefficient blocks are set to zero.
  if (Vlin_star <= tol && Vint_star <= tol) {
    return(list(sL = 0, s = 0, r = NA_real_))
  }

  # If the main-effect target is zero, only the interaction block is scaled.
  if (Vlin_star <= tol) {
    if ((B + G) <= tol) {
      stop("B + G must be > 0 when the interaction target variance is positive.")
    }
    return(list(sL = 0, s = sqrt(Vint_star / (B + G)), r = Inf))
  }

  # If the interaction target is zero, only the main-effect block is scaled.
  if (Vint_star <= tol) {
    if (A <= tol) {
      stop("A must be > 0 when the main-effect target variance is positive.")
    }
    return(list(sL = sqrt(Vlin_star / A), s = 0, r = 0))
  }

  # Without a C/T covariance term, the two blocks can be scaled independently.
  if (abs(D) <= tol) {
    if (A <= tol) {
      stop("A must be > 0.")
    }
    if ((B + G) <= tol) {
      stop("B + G must be > 0.")
    }
    sL <- sqrt(Vlin_star / A)
    s  <- sqrt(Vint_star / (B + G))
    return(list(sL = sL, s = s, r = s / sL))
  }

  # Otherwise, solve the quadratic equation in r = s / sL.
  qa <- Vlin_star * (B + G)
  qb <- D * (Vlin_star - Vint_star)
  qc <- -Vint_star * A

  # A genuinely negative discriminant means that no real ratio r exists.
  disc <- qb^2 - 4 * qa * qc
  if (disc < -tol) {
    stop("No real solution for the ratio r = s / sL.")
  }

  # Tiny negative discriminants can arise from numerical rounding.
  disc <- max(disc, 0)

  # Compute both candidate roots and keep the admissible positive one.
  r1 <- (-qb + sqrt(disc)) / (2 * qa)
  r2 <- (-qb - sqrt(disc)) / (2 * qa)

  admissible <- function(r) is.finite(r) && (r > 0) && ((A + r * D) > tol)
  r <- if (admissible(r1)) r1 else if (admissible(r2)) r2 else NA_real_

  # Stop if neither root yields valid positive scale factors.
  if (!is.finite(r)) {
    stop("No admissible positive solution for r = s / sL.")
  }

  # Once r is known, compute the scale factors.
  sL <- sqrt(Vlin_star / (A + r * D))
  s  <- r * sL

  list(sL = sL, s = s, r = r)
}


# ---- Constant coefficient trajectory -------------------------------------------------------------
# This function repeats Delta_initial for every burn-in and observed wave.

generate_Delta_constant <- function(
  Delta_initial,                       # Initial coefficient matrix Delta_t.
  n_waves,                             # Number of observed waves.
  burn_in = 0L                         # Number of unobserved burn-in waves.
) {

  # Delta_initial must be a matrix because each list element is a Delta_t matrix.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  # Shared trajectory arguments must be valid before the list is built.
  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # The observed trajectory is constant by construction.
  Delta_observed_list <- vector("list", n_waves)
  for (t in seq_len(n_waves)) {
    Delta_observed_list[[t]] <- Delta_initial
  }

  # Burn-in waves are added before the observed trajectory and then labelled.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)
}


# ---- Stepwise coefficient trajectory with a stable rank order ------------------------------------
# This function rescales Delta_initial after step_at without changing coefficient ranks.

generate_Delta_stepwise <- function(
  Delta_initial,                       # Initial coefficient matrix Delta_t.
  n_waves,                             # Number of observed waves.
  burn_in = 0L,                        # Number of unobserved burn-in waves.
  step_at = floor(n_waves / 2) + 1,    # First observed wave receiving the post-step Delta_t.
  R2_old = 0.15,                       # Direct Var(eta_t) before the step.
  R2_new = 0.40                        # Direct Var(eta_t) at and after the step.
) {

  # Delta_initial must be a matrix because it is copied or rescaled over time.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  # Shared trajectory arguments must be valid before the step is placed.
  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # step_at is expressed in observed time, so it must lie between 1 and n_waves.
  if (!is.numeric(step_at) || length(step_at) != 1 || is.na(step_at) ||
      step_at < 1 || step_at > n_waves || step_at != as.integer(step_at)) {
    stop("step_at must be an integer between 1 and n_waves.")
  }

  step_at <- as.integer(step_at)

  # R2_old must be positive because the scale factor divides by it.
  validate_R2(R2_old, name = "R2_old", lower_open = TRUE)

  # R2_new may be zero, which produces a zero post-step coefficient matrix.
  validate_R2(R2_new, name = "R2_new", lower_open = FALSE)

  # Uniformly scaling all coefficients by this factor changes Var(eta_t) by R2_new / R2_old.
  scale_factor <- sqrt(R2_new / R2_old)

  # The observed trajectory changes only at and after step_at.
  Delta_observed_list <- vector("list", n_waves)
  for (t in seq_len(n_waves)) {
    if (t < step_at) {
      Delta_observed_list[[t]] <- Delta_initial
    } else {
      Delta_observed_list[[t]] <- Delta_initial * scale_factor
    }
  }

  # Burn-in waves are added before the observed trajectory and then labelled.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)
}


# ---- Stepwise coefficient trajectory with a mixed post-step pattern ------------------------------
# This function changes coefficient directions after step_at while preserving the old split.

generate_Delta_stepwise_mixture <- function(
  Delta_initial,                       # Initial coefficient matrix Delta_t.
  n_waves,                             # Number of observed waves.
  Omega,                               # Covariance matrix of M = (C', Z', T')'.
  burn_in = 0L,                        # Number of unobserved burn-in waves.
  step_at = floor(n_waves / 2) + 1,    # First observed wave receiving the post-step Delta_t.
  R2_old = 0.15,                       # Direct Var(eta_t) before the step.
  R2_new = 0.40,                       # Direct Var(eta_t) at and after the step.
  lambda_L = 0.50,                     # Mixture weight for the main-effect direction b_L.
  lambda_int = 0.50,                   # Mixture weight for the interaction directions b_2 and b_3.
  seed = NULL,                         # Optional seed for reproducibility.
  tol = 1e-10                          # Numerical tolerance.
) {

  # Shared trajectory arguments must be valid before the post-step matrix is constructed.
  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # step_at is expressed in observed time, not in internal list position.
  if (!is.numeric(step_at) || length(step_at) != 1 || is.na(step_at) ||
      step_at < 1 || step_at > n_waves || step_at != as.integer(step_at)) {
    stop("step_at must be an integer between 1 and n_waves.")
  }

  step_at <- as.integer(step_at)

  # R2_old is the reference direct variance and must be positive.
  validate_R2(R2_old, name = "R2_old", lower_open = TRUE)

  # R2_new may be zero, which produces zero post-step coefficients.
  validate_R2(R2_new, name = "R2_new", lower_open = FALSE)

  # Lambda values determine the degree of rank-order change before rescaling.
  validate_lambda(lambda_L, name = "lambda_L")
  validate_lambda(lambda_int, name = "lambda_int")

  # Omega must be symmetric because it is a covariance matrix.
  if (!is.matrix(Omega) || !isTRUE(all.equal(Omega, t(Omega), tolerance = tol))) {
    stop("Omega must be a symmetric covariance matrix.")
  }

  # The C, Z, and T blocks are identified once and reused for every row.
  blocks <- identify_M_blocks(Delta_initial = Delta_initial, Omega = Omega)

  # The random alternative directions are reproducible when a seed is supplied.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract the covariance blocks needed for the Appendix A equations.
  Omega11 <- Omega[blocks$idx_C, blocks$idx_C, drop = FALSE]
  Omega22 <- if (blocks$has_Z) Omega[blocks$idx_Z, blocks$idx_Z, drop = FALSE] else matrix(0, 0, 0)
  Omega33 <- if (blocks$has_T) Omega[blocks$idx_T, blocks$idx_T, drop = FALSE] else matrix(0, 0, 0)
  Omega13 <- if (blocks$has_T) {
    Omega[blocks$idx_C, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, length(blocks$idx_C), 0)
  }

  # This helper constructs one post-step row from one initial row.
  make_post_row <- function(delta_old_row) {

    # Split the initial row into Appendix A coefficient blocks.
    delta_C_old <- delta_old_row[blocks$idx_C]
    delta_Z_old <- if (blocks$has_Z) delta_old_row[blocks$idx_Z] else numeric(0)
    delta_T_old <- if (blocks$has_T) delta_old_row[blocks$idx_T] else numeric(0)

    # Compute the old effective variance split under Omega.
    old_components <- compute_eta_components(delta_row = delta_old_row, Omega = Omega, blocks = blocks)

    # The supplied initial row must match the stated pre-step direct variance.
    if (abs(old_components$Vtotal - R2_old) > 1e-6) {
      stop("A row of Delta_initial does not match R2_old under the supplied Omega.")
    }

    # The post-step target preserves the old effective main-effect versus interaction split.
    Vlin_star <- (old_components$Vlin / R2_old) * R2_new
    Vint_star <- (old_components$Vint / R2_old) * R2_new

    # Negative effective targets are incompatible with the Appendix A split.
    if (Vlin_star < -tol || Vint_star < -tol) {
      stop("The implied Appendix A variance split contains a negative component.")
    }

    Vlin_star <- max(Vlin_star, 0)
    Vint_star <- max(Vint_star, 0)

    # Sample random alternative directions for the same C, Z, and T blocks.
    b_L_new <- rnorm(length(delta_C_old))
    b_2_new <- if (blocks$has_Z) rnorm(length(delta_Z_old)) else numeric(0)
    b_3_new <- if (blocks$has_T) rnorm(length(delta_T_old)) else numeric(0)

    # The mixed vectors are directions only and are rescaled below.
    b_L_mix <- (1 - lambda_L) * delta_C_old + lambda_L * b_L_new
    b_2_mix <- if (blocks$has_Z) {
      (1 - lambda_int) * delta_Z_old + lambda_int * b_2_new
    } else {
      numeric(0)
    }
    b_3_mix <- if (blocks$has_T) {
      (1 - lambda_int) * delta_T_old + lambda_int * b_3_new
    } else {
      numeric(0)
    }

    # Compute the unscaled Appendix A quantities for the mixed direction.
    A <- qform(b_L_mix, Omega11)
    B <- if (blocks$has_Z) qform(b_2_mix, Omega22) else 0
    G <- if (blocks$has_T) qform(b_3_mix, Omega33) else 0
    D <- if (blocks$has_T) as.numeric(t(b_L_mix) %*% Omega13 %*% b_3_mix) else 0

    # Solve for the scale factors sL and s.
    scales <- solve_appendix_scales(
      A = A,
      B = B,
      G = G,
      D = D,
      Vlin_star = Vlin_star,
      Vint_star = Vint_star,
      tol = tol
    )

    # Apply the scale factors to obtain the post-step coefficient blocks.
    delta_C_post <- scales$sL * b_L_mix
    delta_Z_post <- if (blocks$has_Z) scales$s * b_2_mix else numeric(0)
    delta_T_post <- if (blocks$has_T) scales$s * b_3_mix else numeric(0)

    # Reconstruct the row in the same order as M = (C', Z', T')'.
    delta_post_row <- numeric(length(delta_old_row))
    delta_post_row[blocks$idx_C] <- delta_C_post
    if (blocks$has_Z) delta_post_row[blocks$idx_Z] <- delta_Z_post
    if (blocks$has_T) delta_post_row[blocks$idx_T] <- delta_T_post

    delta_post_row
  }

  # Build the post-step coefficient matrix row by row.
  Delta_post <- Delta_initial
  for (row in seq_len(nrow(Delta_initial))) {
    Delta_post[row, ] <- make_post_row(Delta_initial[row, ])
  }

  rownames(Delta_post) <- rownames(Delta_initial)
  colnames(Delta_post) <- colnames(Delta_initial)

  # The observed trajectory changes only at and after step_at.
  Delta_observed_list <- vector("list", n_waves)
  for (t in seq_len(n_waves)) {
    if (t < step_at) {
      Delta_observed_list[[t]] <- Delta_initial
    } else {
      Delta_observed_list[[t]] <- Delta_post
    }
  }

  # Burn-in waves are added before the observed trajectory and then labelled.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  Delta_list <- name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)

  # The post-step matrix is stored for inspection without changing the list structure.
  attr(Delta_list, "Delta_post") <- Delta_post

  Delta_list
}
