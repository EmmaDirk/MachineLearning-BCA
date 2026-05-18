# ------------------------------------------------------------------------------------------------
#
# This script contains helper functions for constructing trajectories of coefficient matrices
# Delta_t. Each Delta_t determines how the expanded baseline-confounder vector M contributes
# directly to X_t and Y_t in the data-generating mechanism.
#
# The notation follows Appendix A:
#
#   eta_t = Delta_t M,
#   M     = (C', Z', T')'.
#
# C contains the base confounders, Z contains standardized two-way interaction terms, and T
# contains standardized three-way interaction terms. The covariance matrix of M is denoted Omega.
#
# Three trajectory types are implemented:
#
# 1. constant effects
#    Delta_t is kept fixed for all generated waves.
#
# 2. stepwise effects
#    Delta_t is kept fixed before a chosen observed wave and then multiplied by one common scale
#    factor. This changes Var(eta_t) from R2_old to R2_new while preserving the coefficient pattern.
#
# 3. stepwise mixture effects
#    The post-step Delta_t is constructed as a mixture of the initial coefficient pattern and a
#    randomly sampled alternative pattern. The mixed direction is then rescaled using the same
#    Appendix A logic as the initial coefficient sampler. This allows the rank order of the
#    coefficients to change while still targeting the requested direct variance of eta_t.
#
# Important caution:
#
# R2_old and R2_new refer to the direct variance of eta_t at the wave where the corresponding
# Delta_t is applied. They do not necessarily equal the total baseline-confounder-related variance
# in X_t or Y_t after autoregressive and cross-lagged propagation. If the data-generating mechanism
# includes a burn-in period and the retained waves are interpreted after the process has moved
# toward stationarity, the realised confounder-related variance in those retained waves can differ
# from R2_old and R2_new. These functions only construct Delta_t trajectories; they do not verify
# the propagated variance of the full longitudinal system.
#
# burn_in is interpreted as the number of unobserved waves to prepend before the observed waves.
# step_at is interpreted in observed time. For example, if step_at = 3 and burn_in = 20, then the
# post-step coefficient matrix starts at position 23 in the returned list.
#
# ------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------
# Shared input validation
# ------------------------------------------------------------------------------------------------

validate_Delta_trajectory_inputs <- function(n_waves, burn_in) {

  # n_waves is the number of observed waves requested by the user.
  if (!is.numeric(n_waves) || length(n_waves) != 1 ||
      n_waves < 1 || n_waves != as.integer(n_waves)) {
    stop("n_waves must be a positive integer.")
  }

  # burn_in is the number of unobserved waves prepended to the observed trajectory.
  if (!is.numeric(burn_in) || length(burn_in) != 1 ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("burn_in must be a non-negative integer.")
  }
}


validate_R2 <- function(R2, name, lower_open = FALSE) {

  # R2-type arguments are interpreted as variances of eta_t when Var(X_t) or Var(Y_t) is 1.
  if (!is.numeric(R2) || length(R2) != 1 || is.na(R2)) {
    stop(paste0(name, " must be a single numeric value."))
  }

  # Some arguments, such as R2_old in a rescaling ratio, must be strictly positive.
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

  # Lambda controls how strongly the post-step direction is pulled toward a random alternative.
  if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) ||
      lambda < 0 || lambda > 1) {
    stop(paste0(name, " must be a single number in [0, 1]."))
  }
}


# ------------------------------------------------------------------------------------------------
# Helper for adding burn-in waves
# ------------------------------------------------------------------------------------------------

prepend_burn_in <- function(Delta_observed_list, Delta_initial, burn_in) {

  # If there is no burn-in, the observed trajectory is already the full returned trajectory.
  if (burn_in == 0L) {
    return(Delta_observed_list)
  }

  # During burn-in, use copies of the initial coefficient matrix. This means that any stepwise
  # change is placed only in the observed part of the trajectory.
  Delta_burn_in_list <- vector("list", burn_in)

  # Fill every burn-in position with the initial Delta_t.
  for (b in seq_len(burn_in)) {
    Delta_burn_in_list[[b]] <- Delta_initial
  }

  # Return the internal trajectory: first burn-in waves, then observed waves.
  c(Delta_burn_in_list, Delta_observed_list)
}


name_Delta_trajectory <- function(Delta_list, burn_in) {

  # Names are internal positions in the returned list. They are not meant to encode any special
  # substantive time labels such as -20, ..., 0, 1, ..., 5.
  names(Delta_list) <- paste0("t", seq_along(Delta_list))

  # Mark which positions are burn-in and which are observed. This is useful when burn_in > 0.
  attr(Delta_list, "burn_in") <- burn_in
  attr(Delta_list, "n_total") <- length(Delta_list)
  attr(Delta_list, "n_observed") <- length(Delta_list) - burn_in

  Delta_list
}


# ------------------------------------------------------------------------------------------------
# Helper for identifying the blocks of M = (C', Z', T')'
# ------------------------------------------------------------------------------------------------

identify_M_blocks <- function(Delta_initial, Omega) {

  # Delta_initial must be a matrix with one row per outcome and one column per feature in M.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  # Omega must be the covariance matrix of M, in the same feature order as the columns of
  # Delta_initial.
  if (!is.matrix(Omega)) {
    stop("Omega must be a matrix.")
  }

  # Delta_initial and Omega must refer to the same number of features.
  if (ncol(Delta_initial) != nrow(Omega) || nrow(Omega) != ncol(Omega)) {
    stop("Delta_initial and Omega are not conformable.")
  }

  # Feature names are required because the blocks C, Z, and T are identified from their prefixes.
  if (is.null(colnames(Delta_initial))) {
    stop("Delta_initial must have column names following the C/Z/T naming scheme.")
  }

  if (is.null(rownames(Omega)) || is.null(colnames(Omega))) {
    stop("Omega must have row and column names following the C/Z/T naming scheme.")
  }

  # The names and order of Omega must match Delta_initial exactly.
  if (!identical(colnames(Delta_initial), rownames(Omega)) ||
      !identical(colnames(Delta_initial), colnames(Omega))) {
    stop("The feature names of Delta_initial and Omega must match exactly and be in the same order.")
  }

  # C names identify the base-confounder block.
  idx_C <- grep("^C[0-9]+$", colnames(Delta_initial))

  # Z names identify standardized two-way interaction terms.
  idx_Z <- grep("^Z[0-9]+:[0-9]+$", colnames(Delta_initial))

  # T names identify standardized three-way interaction terms.
  idx_T <- grep("^T[0-9]+:[0-9]+:[0-9]+$", colnames(Delta_initial))

  # There must be at least one base confounder.
  if (length(idx_C) == 0L) {
    stop("No C block found. Column names should include C1, C2, ..., Ck.")
  }

  # The expected Appendix A order is C first, then Z, then T.
  idx_expected <- c(idx_C, idx_Z, idx_T)

  if (!identical(idx_expected, seq_along(idx_expected))) {
    stop("Columns must be ordered as M = (C', Z', T')'.")
  }

  # Return block indices and logical flags used by the rescaling functions.
  list(
    idx_C = idx_C,
    idx_Z = idx_Z,
    idx_T = idx_T,
    has_Z = length(idx_Z) > 0L,
    has_T = length(idx_T) > 0L,
    has_int = (length(idx_Z) + length(idx_T)) > 0L
  )
}


# ------------------------------------------------------------------------------------------------
# Quadratic forms and Appendix A variance decomposition
# ------------------------------------------------------------------------------------------------

qform <- function(x, Omega_block) {

  # Compute x' Omega_block x and return it as a scalar.
  as.numeric(t(x) %*% Omega_block %*% x)
}


compute_eta_components <- function(delta_row, Omega, blocks) {

  # Split one row of Delta_t into the coefficient blocks corresponding to C, Z, and T.
  delta_C <- delta_row[blocks$idx_C]
  delta_Z <- if (blocks$has_Z) delta_row[blocks$idx_Z] else numeric(0)
  delta_T <- if (blocks$has_T) delta_row[blocks$idx_T] else numeric(0)

  # Extract the covariance blocks from Omega.
  Omega11 <- Omega[blocks$idx_C, blocks$idx_C, drop = FALSE]
  Omega22 <- if (blocks$has_Z) {
    Omega[blocks$idx_Z, blocks$idx_Z, drop = FALSE]
  } else {
    matrix(0, 0, 0)
  }
  Omega33 <- if (blocks$has_T) {
    Omega[blocks$idx_T, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, 0, 0)
  }
  Omega13 <- if (blocks$has_T) {
    Omega[blocks$idx_C, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, length(blocks$idx_C), 0)
  }

  # A is the variance contribution of the C block before the C/T cross term is added.
  A <- qform(delta_C, Omega11)

  # B is the variance contribution of the Z block.
  B <- if (blocks$has_Z) qform(delta_Z, Omega22) else 0

  # G is the variance contribution of the T block.
  G <- if (blocks$has_T) qform(delta_T, Omega33) else 0

  # D is the covariance term between C and T.
  D <- if (blocks$has_T) {
    as.numeric(t(delta_C) %*% Omega13 %*% delta_T)
  } else {
    0
  }

  # Appendix A splits the total C/T cross term symmetrically between the main-effect and
  # interaction components.
  Vlin <- A + D
  Vint <- B + G + D

  # The total direct variance of eta_t is the sum of both effective components.
  Vtotal <- Vlin + Vint

  list(
    A = A,
    B = B,
    G = G,
    D = D,
    Vlin = Vlin,
    Vint = Vint,
    Vtotal = Vtotal
  )
}


# ------------------------------------------------------------------------------------------------
# Solve Appendix A scale factors
# ------------------------------------------------------------------------------------------------

solve_appendix_scales <- function(A, B, G, D, Vlin_star, Vint_star, tol = 1e-10) {

  # This helper solves for the same two scale factors used in Appendix A:
  #
  # - sL scales the main-effect direction b_L;
  # - s  scales the interaction directions b_2 and b_3.
  #
  # The unscaled quantities are:
  #
  # A = b_L' Omega11 b_L,
  # B = b_2' Omega22 b_2,
  # G = b_3' Omega33 b_3,
  # D = b_L' Omega13 b_3.
  #
  # The target effective components are:
  #
  # Vlin_star = sL^2 A       + sL s D,
  # Vint_star = s^2 (B + G) + sL s D.

  # If the target total is zero, both scale factors are zero.
  if (Vlin_star <= tol && Vint_star <= tol) {
    return(list(sL = 0, s = 0, r = NA_real_))
  }

  # If the target main-effect component is zero, only the interaction block remains.
  if (Vlin_star <= tol) {

    if ((B + G) <= tol) {
      stop("B + G must be > 0 when the interaction target variance is positive.")
    }

    return(list(sL = 0, s = sqrt(Vint_star / (B + G)), r = Inf))
  }

  # If the target interaction component is zero, only the main-effect block remains.
  if (Vint_star <= tol) {

    if (A <= tol) {
      stop("A must be > 0 when the main-effect target variance is positive.")
    }

    return(list(sL = sqrt(Vlin_star / A), s = 0, r = 0))
  }

  # If there is no C/T covariance term, the two blocks can be scaled independently.
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

  # Otherwise solve the quadratic equation in r = s / sL:
  #
  # Vlin_star (B + G) r^2 + D (Vlin_star - Vint_star) r - Vint_star A = 0.
  qa <- Vlin_star * (B + G)
  qb <- D * (Vlin_star - Vint_star)
  qc <- -Vint_star * A

  # Compute the discriminant.
  disc <- qb^2 - 4 * qa * qc

  # A genuinely negative discriminant means that no real ratio r exists.
  if (disc < -tol) {
    stop("No real solution for the ratio r = s / sL.")
  }

  # Guard against tiny negative values caused by numerical rounding.
  disc <- max(disc, 0)

  # Compute both candidate roots.
  r1 <- (-qb + sqrt(disc)) / (2 * qa)
  r2 <- (-qb - sqrt(disc)) / (2 * qa)

  # A root is admissible if it is positive and makes the denominator for sL positive.
  admissible <- function(r) is.finite(r) && (r > 0) && ((A + r * D) > tol)

  # Select the admissible positive root.
  r <- if (admissible(r1)) r1 else if (admissible(r2)) r2 else NA_real_

  # Stop if neither root is admissible.
  if (!is.finite(r)) {
    stop("No admissible positive solution for r = s / sL.")
  }

  # Once r is known, compute the scale factors.
  sL <- sqrt(Vlin_star / (A + r * D))
  s  <- r * sL

  list(sL = sL, s = s, r = r)
}


# ------------------------------------------------------------------------------------------------
# 1) Constant coefficient trajectory
# ------------------------------------------------------------------------------------------------

generate_Delta_constant <- function(
  Delta_initial,                       # initial coefficient matrix Delta_t
  n_waves,                             # number of observed waves
  burn_in = 0L                         # number of unobserved burn-in waves
) {

  # Check shared inputs.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  # Convert integer-like inputs to integer storage.
  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # Create the observed part of the trajectory.
  Delta_observed_list <- vector("list", n_waves)

  # Keep Delta_t fixed across all observed waves.
  for (t in seq_len(n_waves)) {
    Delta_observed_list[[t]] <- Delta_initial
  }

  # Add the burn-in part before the observed part.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  # Name the internal positions and attach simple metadata.
  Delta_list <- name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)

  return(Delta_list)
}


# ------------------------------------------------------------------------------------------------
# 2) Stepwise coefficient trajectory with preserved coefficient pattern
# ------------------------------------------------------------------------------------------------

generate_Delta_stepwise <- function(
  Delta_initial,                       # initial coefficient matrix Delta_t
  n_waves,                             # number of observed waves
  burn_in = 0L,                        # number of unobserved burn-in waves
  step_at = floor(n_waves / 2) + 1,    # first observed wave receiving the post-step Delta_t
  R2_old = 0.15,                       # direct Var(eta_t) before the step
  R2_new = 0.40                        # direct Var(eta_t) at and after the step
) {

  # Check shared inputs.
  if (!is.matrix(Delta_initial)) {
    stop("Delta_initial must be a matrix.")
  }

  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  # Convert integer-like inputs to integer storage.
  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # The step is specified in observed time, not internal list position.
  if (!is.numeric(step_at) || length(step_at) != 1 ||
      step_at < 1 || step_at > n_waves || step_at != as.integer(step_at)) {
    stop("step_at must be an integer between 1 and n_waves.")
  }

  step_at <- as.integer(step_at)

  # R2_old must be positive because the scale factor divides by it.
  validate_R2(R2_old, name = "R2_old", lower_open = TRUE)

  # R2_new may be zero, in which case post-step coefficients are set to zero.
  validate_R2(R2_new, name = "R2_new", lower_open = FALSE)

  # Multiplying every coefficient by sqrt(R2_new / R2_old) multiplies Var(eta_t)
  # by R2_new / R2_old. This preserves coefficient signs, rank order, and the relative
  # Appendix A variance split.
  scale_factor <- sqrt(R2_new / R2_old)

  # Create the observed part of the trajectory.
  Delta_observed_list <- vector("list", n_waves)

  # Fill the observed trajectory.
  for (t in seq_len(n_waves)) {

    # Before the step, use the initial coefficient matrix.
    if (t < step_at) {
      Delta_observed_list[[t]] <- Delta_initial

    # At and after the step, use the uniformly rescaled coefficient matrix.
    } else {
      Delta_observed_list[[t]] <- Delta_initial * scale_factor
    }
  }

  # Add the burn-in part before the observed part.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  # Name the internal positions and attach simple metadata.
  Delta_list <- name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)

  return(Delta_list)
}


# ------------------------------------------------------------------------------------------------
# 3) Stepwise coefficient trajectory with a mixed post-step pattern
# ------------------------------------------------------------------------------------------------

generate_Delta_stepwise_mixture <- function(
  Delta_initial,                       # initial coefficient matrix Delta_t
  n_waves,                             # number of observed waves
  Omega,                               # covariance matrix of M = (C', Z', T')'
  burn_in = 0L,                        # number of unobserved burn-in waves
  step_at = floor(n_waves / 2) + 1,    # first observed wave receiving the post-step Delta_t
  R2_old = 0.15,                       # direct Var(eta_t) before the step
  R2_new = 0.40,                       # direct Var(eta_t) at and after the step
  lambda_L = 0.50,                     # mixture weight for the main-effect direction b_L
  lambda_int = 0.50,                   # mixture weight for interaction directions b_2 and b_3
  seed = NULL,                         # optional seed for reproducibility
  tol = 1e-10                          # numerical tolerance
) {

  # Check shared inputs.
  validate_Delta_trajectory_inputs(n_waves = n_waves, burn_in = burn_in)

  # Convert integer-like inputs to integer storage.
  n_waves <- as.integer(n_waves)
  burn_in <- as.integer(burn_in)

  # The step is specified in observed time, not internal list position.
  if (!is.numeric(step_at) || length(step_at) != 1 ||
      step_at < 1 || step_at > n_waves || step_at != as.integer(step_at)) {
    stop("step_at must be an integer between 1 and n_waves.")
  }

  step_at <- as.integer(step_at)

  # R2_old must be positive because the old direct variance is used as a reference.
  validate_R2(R2_old, name = "R2_old", lower_open = TRUE)

  # R2_new may be zero, in which case post-step coefficients are set to zero.
  validate_R2(R2_new, name = "R2_new", lower_open = FALSE)

  # Lambda values control how much the direction changes before rescaling.
  validate_lambda(lambda_L, name = "lambda_L")
  validate_lambda(lambda_int, name = "lambda_int")

  # Omega must be symmetric because it is a covariance matrix.
  if (!is.matrix(Omega) || !isTRUE(all.equal(Omega, t(Omega), tolerance = tol))) {
    stop("Omega must be a symmetric covariance matrix.")
  }

  # Identify the C, Z, and T blocks in the columns of Delta_initial and in Omega.
  blocks <- identify_M_blocks(Delta_initial = Delta_initial, Omega = Omega)

  # Optional reproducibility for the random alternative directions.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract covariance blocks once, because they are reused for every row.
  Omega11 <- Omega[blocks$idx_C, blocks$idx_C, drop = FALSE]
  Omega22 <- if (blocks$has_Z) {
    Omega[blocks$idx_Z, blocks$idx_Z, drop = FALSE]
  } else {
    matrix(0, 0, 0)
  }
  Omega33 <- if (blocks$has_T) {
    Omega[blocks$idx_T, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, 0, 0)
  }
  Omega13 <- if (blocks$has_T) {
    Omega[blocks$idx_C, blocks$idx_T, drop = FALSE]
  } else {
    matrix(0, length(blocks$idx_C), 0)
  }

  # This helper constructs one post-step coefficient row from one initial row.
  make_post_row <- function(delta_old_row) {

    # Split the initial coefficient row into the Appendix A blocks.
    delta_C_old <- delta_old_row[blocks$idx_C]
    delta_Z_old <- if (blocks$has_Z) delta_old_row[blocks$idx_Z] else numeric(0)
    delta_T_old <- if (blocks$has_T) delta_old_row[blocks$idx_T] else numeric(0)

    # Compute the Appendix A decomposition of the old direct variance.
    old_components <- compute_eta_components(
      delta_row = delta_old_row,
      Omega = Omega,
      blocks = blocks
    )

    # Check that the supplied initial row actually targets R2_old under Omega.
    if (abs(old_components$Vtotal - R2_old) > 1e-6) {
      stop("A row of Delta_initial does not match R2_old under the supplied Omega.")
    }

    # Preserve the old effective Appendix A variance split for this row.
    # This keeps the same relative main-effect and interaction contributions while
    # changing the total direct variance from R2_old to R2_new.
    Vlin_star <- (old_components$Vlin / R2_old) * R2_new
    Vint_star <- (old_components$Vint / R2_old) * R2_new

    # The effective components should be non-negative. Tiny negative values can arise from
    # numerical rounding, but substantively negative targets are incompatible with the split.
    if (Vlin_star < -tol || Vint_star < -tol) {
      stop("The implied Appendix A variance split contains a negative component.")
    }

    Vlin_star <- max(Vlin_star, 0)
    Vint_star <- max(Vint_star, 0)

    # Sample random alternative directions for the same C, Z, and T blocks.
    b_L_new <- rnorm(length(delta_C_old))
    b_2_new <- if (blocks$has_Z) rnorm(length(delta_Z_old)) else numeric(0)
    b_3_new <- if (blocks$has_T) rnorm(length(delta_T_old)) else numeric(0)

    # Mix the initial direction with the random alternative direction. The mixed vector is only
    # a direction; it is rescaled below to hit the requested variance targets.
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
    D <- if (blocks$has_T) {
      as.numeric(t(b_L_mix) %*% Omega13 %*% b_3_mix)
    } else {
      0
    }

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

    # Apply the Appendix A scale factors to obtain the post-step blocks.
    delta_C_post <- scales$sL * b_L_mix
    delta_Z_post <- if (blocks$has_Z) scales$s * b_2_mix else numeric(0)
    delta_T_post <- if (blocks$has_T) scales$s * b_3_mix else numeric(0)

    # Reconstruct the full row in the same order as M = (C', Z', T')'.
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

  # Preserve row and column names.
  rownames(Delta_post) <- rownames(Delta_initial)
  colnames(Delta_post) <- colnames(Delta_initial)

  # Create the observed part of the trajectory.
  Delta_observed_list <- vector("list", n_waves)

  # Fill the observed trajectory.
  for (t in seq_len(n_waves)) {

    # Before the step, use the initial coefficient matrix.
    if (t < step_at) {
      Delta_observed_list[[t]] <- Delta_initial

    # At and after the step, use the mixed-and-rescaled coefficient matrix.
    } else {
      Delta_observed_list[[t]] <- Delta_post
    }
  }

  # Add the burn-in part before the observed part.
  Delta_list <- prepend_burn_in(
    Delta_observed_list = Delta_observed_list,
    Delta_initial = Delta_initial,
    burn_in = burn_in
  )

  # Name the internal positions and attach simple metadata.
  Delta_list <- name_Delta_trajectory(Delta_list = Delta_list, burn_in = burn_in)

  # Store the post-step matrix as an attribute for direct inspection.
  attr(Delta_list, "Delta_post") <- Delta_post

  return(Delta_list)
}
