# =================================================================================================
#
# This script defines a function for sampling the coefficient matrix Delta_t. Delta_t
# determines how the baseline-confounder term vector M contributes directly to X_t
# and Y_t in the data-generating mechanism.
#
# For one generic outcome W_t, the direct baseline-confounder contribution is
#
#   eta_t = delta_t' M,
#
# where M contains the baseline confounders and, optionally, standardized two-way and
# three-way interaction terms. The function samples delta_t so that
#
#   Var(eta_t) = R2_total.
#
# If W_t is generated with variance 1, this variance can be interpreted as the direct
# R^2 of the baseline-confounder term at the wave where the sampled coefficients are
# applied.
#
# Caution: this function scales only the direct contribution eta_t. In longitudinal
# simulations with autoregressive and cross-lagged paths, baseline-confounder effects
# can propagate across waves. If a burn-in period is simulated and discarded, R2_total
# should therefore not be read as the total baseline-confounder-related variance in the
# retained measurement waves. It is the variance of eta_t at the scaling point, not the
# total propagated confounder effect after the system has evolved.
#
# The argument k gives the number of base confounders:
#
#   C1, C2, ..., Ck.
#
# The expanded vector M has up to three blocks:
#
#   M = (C', Z', T')',
#
# where C contains the base confounders, Z contains standardized two-way interactions,
# and T contains standardized three-way interactions.
#
# Interaction terms are included only when requested:
#
# - include_2way adds all possible two-way interactions, such as Z1:2, Z1:3, ...
# - include_3way adds all possible three-way interactions, such as T1:2:4, ...
#
# The function first samples random coefficient directions from a normal distribution.
# These directions determine the relative signs and magnitudes of the coefficients.
# It then rescales the directions so that the main-effect and interaction blocks explain
# the requested shares of R2_total.
#
# When three-way interactions are included, the main-effect block C and the three-way
# block T can be correlated. Following Appendix A, the resulting cross term is split
# symmetrically: one half is assigned to the main-effect contribution and one half to
# the interaction contribution.
#
# The procedure relies on the following assumptions:
#
# - the base confounders are centered, so E(C_i) = 0;
# - the base confounders are standardized, so Var(C_i) = 1;
# - the base confounders are jointly Gaussian;
# - Omega11 is the covariance matrix of the base confounders;
# - if R2_total is interpreted as an R^2, W_t is generated with variance 1.
#
# The interaction terms themselves are standardized products of Gaussian variables.
# They are centered and have variance 1, but they are not normally distributed.
#
# =================================================================================================

sample_delta_t <- function(
  k,                             # number of base confounders in C
  Omega11,                       # k x k covariance matrix of C
  R2_total,                      # target Var(eta_t); equals direct R^2 if Var(W_t) = 1
  rho_int = 0,                   # share of R2_total assigned to interaction terms
  include_2way = FALSE,          # include all possible standardized two-way terms Z
  include_3way = FALSE,          # include all possible standardized three-way terms T
  force_positive = FALSE,        # sample only non-negative coefficient directions
  min_abs = 0,                   # minimum allowed absolute final coefficient size
  max_abs = 1,                   # maximum allowed absolute final coefficient size
  max_tries = 100000,            # maximum number of resampling attempts
  eig_tol = 1e-10                # tolerance for the positive-semidefinite check
) {

  # ---- input checks ------------------------------------------------------------------------------
  #
  # This section checks whether the input arguments are compatible with the
  # coefficient-sampling procedure described in Appendix A.

  # Omega11 must be a k x k matrix because it represents the covariance matrix
  # of the k base confounders in C.
  if (!is.matrix(Omega11) || nrow(Omega11) != k || ncol(Omega11) != k)
    stop("Omega11 must be a k x k matrix.")

  # Omega11 must be symmetric because covariance matrices are symmetric.
  if (!isTRUE(all.equal(Omega11, t(Omega11), tolerance = 1e-12)))
    stop("Omega11 must be symmetric.")

  # Omega11 must be positive semidefinite because all variances implied by a
  # covariance matrix must be non-negative.
  ev <- eigen(Omega11, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < -eig_tol))
    stop("Omega11 must be positive semidefinite.")

  # The diagonal of Omega11 must be 1 because the base confounders are assumed
  # to be standardized.
  if (!isTRUE(all.equal(diag(Omega11), rep(1, k), tolerance = 1e-12)))
    stop("Omega11 must have 1 on the diagonal because C is standardized.")

  # R2_total is interpreted as a variance share and must therefore lie between
  # 0 and 1.
  if (R2_total < 0 || R2_total > 1)
    stop("R2_total must be between 0 and 1.")

  # rho_int is the target share of R2_total assigned to the interaction block.
  # It must also lie between 0 and 1.
  if (rho_int < 0 || rho_int > 1)
    stop("rho_int must be between 0 and 1.")

  # Two-way interactions require at least two base confounders.
  if (include_2way && k < 2)
    stop("include_2way requires k >= 2.")

  # Three-way interactions require at least three base confounders.
  if (include_3way && k < 3)
    stop("include_3way requires k >= 3.")

  # If interaction terms are included, their target variance share must be
  # positive. Otherwise the interaction coefficients would be scaled to zero.
  if ((include_2way || include_3way) && rho_int == 0)
    stop("rho_int must be > 0 when interaction terms are included.")

  # If rho_int is positive, at least one interaction block must be included.
  if (rho_int > 0 && !(include_2way || include_3way))
    stop("rho_int is positive but no interaction terms are included.")

  # The lower bound for absolute coefficient sizes cannot be negative.
  if (min_abs < 0)
    stop("min_abs must be >= 0.")

  # The upper bound must be at least as large as the lower bound.
  if (max_abs < min_abs)
    stop("max_abs must be >= min_abs.")

  # ---- dimensions --------------------------------------------------------------------------------
  #
  # This section determines the number of two-way and three-way interaction
  # terms that appear in the expanded vector M.

  # The number of possible two-way interaction terms is k choose 2.
  m2 <- if (include_2way) choose(k, 2) else 0

  # The number of possible three-way interaction terms is k choose 3.
  m3 <- if (include_3way) choose(k, 3) else 0

  # The total number of interaction terms is the length of Z plus the length of T.
  m_int <- m2 + m3

  # The target main-effect contribution is (1 - rho_int) times R2_total.
  Vlin_star <- (1 - rho_int) * R2_total

  # The target interaction contribution is rho_int times R2_total.
  VNL_star <- rho_int * R2_total

  # ---- helper functions --------------------------------------------------------------------------
  #
  # These functions construct the index sets needed for the interaction terms.

  # Return all index pairs used to form Z.
  # For k = 3, this returns rows (1, 2), (1, 3), and (2, 3).
  all_pairs <- function(k) t(combn(seq_len(k), 2))

  # Return all index triples used to form T.
  # For k = 4, this includes rows such as (1, 2, 3) and (1, 2, 4).
  all_triples <- function(k) t(combn(seq_len(k), 3))

  # ---- covariance objects ------------------------------------------------------------------------
  #
  # This section computes the covariance blocks of the expanded vector M.
  #
  # The blocks follow the notation in Appendix A:
  #
  # - Omega11 is Cov(C, C), the covariance matrix of the base confounders;
  # - Omega22 is Cov(Z, Z), the covariance matrix of two-way interactions;
  # - Omega33 is Cov(T, T), the covariance matrix of three-way interactions;
  # - Omega13 is Cov(C, T), the covariance matrix between C and T.
  #
  # Because the base confounders are centered and jointly Gaussian, odd-order
  # moments are zero. Therefore:
  #
  # - Omega12 = Cov(C, Z) = 0;
  # - Omega23 = Cov(Z, T) = 0.

  # Helper:
  # For centered jointly Gaussian variables, Isserlis' theorem says that a
  # sixth-order moment is the sum over all 15 pairings of the variables. The
  # pairings are precomputed once and reused below.
  pairings6 <- local({

    # Recursively construct all pairings of a vector of positions.
    rec <- function(v) {

      # If no positions remain, the current pairing is complete.
      if (length(v) == 0) return(list(list()))

      # If two positions remain, they must be paired with each other.
      if (length(v) == 2) return(list(list(c(v[1], v[2]))))

      # Otherwise, take the first position and pair it with each possible
      # remaining position.
      first <- v[1]
      out <- list()

      for (m in 2:length(v)) {

        # Remove the two positions that have just been paired.
        rest <- v[-c(1, m)]

        # Recursively pair all remaining positions.
        sub <- rec(rest)

        # Add the newly created pair to each recursive pairing.
        for (s in sub) {
          out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
        }
      }

      out
    }

    # There are 15 pairings of six positions.
    rec(1:6)
  })

  # Helper:
  # Compute E[X1 X2 X3 X4 X5 X6] for centered jointly Gaussian variables.
  # The argument idx contains the indices of the base confounders appearing in
  # the product, and Omega11 supplies their pairwise covariances.
  sixth_moment_gaussian <- function(idx, Omega11) {

    # For each of the 15 pairings, multiply the three covariance terms. Then
    # sum these products over pairings.
    sum(vapply(
      pairings6,
      function(pr) {
        prod(vapply(
          pr,
          function(pa) Omega11[idx[pa[1]], idx[pa[2]]],
          numeric(1)
        ))
      },
      numeric(1)
    ))
  }

  # Analytic covariance between two standardized two-way interactions.
  # This computes Cov(Z_ij, Z_ab), corresponding to Equation A49.
  cov_Z2 <- function(i, j, a, b, Omega11) {

    # The numerator is the covariance between the centered raw products
    # C_i C_j and C_a C_b.
    num <- Omega11[i, a] * Omega11[j, b] +
           Omega11[i, b] * Omega11[j, a]

    # The denominator standardizes both two-way products to variance 1.
    den <- sqrt((1 + Omega11[i, j]^2) *
                (1 + Omega11[a, b]^2))

    num / den
  }

  # Analytic covariance between two standardized three-way interactions.
  # This computes Cov(T_ijl, T_abc), corresponding to Equation A51.
  cov_T3 <- function(i, j, l, a, b, c, Omega11) {

    # The numerator is E[C_i C_j C_l C_a C_b C_c]. It is evaluated using
    # Isserlis' theorem through sixth_moment_gaussian().
    num <- sixth_moment_gaussian(c(i, j, l, a, b, c), Omega11)

    # The first denominator term is v_ijl, the variance of C_i C_j C_l.
    v_ijl <- 1 +
      2 * Omega11[i, j]^2 +
      2 * Omega11[i, l]^2 +
      2 * Omega11[j, l]^2 +
      8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]

    # The second denominator term is v_abc, the variance of C_a C_b C_c.
    v_abc <- 1 +
      2 * Omega11[a, b]^2 +
      2 * Omega11[a, c]^2 +
      2 * Omega11[b, c]^2 +
      8 * Omega11[a, b] * Omega11[a, c] * Omega11[b, c]

    # Dividing by sqrt(v_ijl * v_abc) standardizes both three-way products.
    num / sqrt(v_ijl * v_abc)
  }

  # Analytic covariance between a main effect and a standardized three-way
  # interaction. This computes Cov(C_a, T_ijl), corresponding to Equation A50.
  cov_C_T <- function(a, i, j, l, Omega11) {

    # The numerator is E[C_a C_i C_j C_l], expanded by Isserlis' theorem.
    num <- Omega11[a, i] * Omega11[j, l] +
           Omega11[a, j] * Omega11[i, l] +
           Omega11[a, l] * Omega11[i, j]

    # The denominator is sqrt(v_ijl), the standard deviation of C_i C_j C_l.
    den <- sqrt(
      1 +
      2 * Omega11[i, j]^2 +
      2 * Omega11[i, l]^2 +
      2 * Omega11[j, l]^2 +
      8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]
    )

    num / den
  }

  # Build Omega22, the covariance matrix of the standardized two-way
  # interaction block Z.
  build_Omega22 <- function(Omega11) {

    # List all pairs that define entries of Z.
    P <- all_pairs(nrow(Omega11))

    # Number of two-way interaction terms.
    m2 <- nrow(P)

    # Initialize the covariance matrix of Z.
    Omega22 <- matrix(0, m2, m2)

    # If there are no pairs, return the empty matrix.
    if (m2 == 0) return(Omega22)

    # Each entry of Z is standardized, so the diagonal is 1.
    diag(Omega22) <- 1

    # Fill the off-diagonal covariances analytically.
    for (row in seq_len(m2)) {

      # Indices of the first two-way interaction.
      i <- P[row, 1]
      j <- P[row, 2]

      for (col in seq_len(m2)) {

        # The diagonal has already been set to 1.
        if (row == col) next

        # Indices of the second two-way interaction.
        a <- P[col, 1]
        b <- P[col, 2]

        # Fill Cov(Z_ij, Z_ab).
        Omega22[row, col] <- cov_Z2(i, j, a, b, Omega11)
      }
    }

    Omega22
  }

  # Build Omega33, the covariance matrix of the standardized three-way
  # interaction block T.
  build_Omega33 <- function(Omega11) {

    # List all triples that define entries of T.
    Tr <- all_triples(nrow(Omega11))

    # Number of three-way interaction terms.
    m3 <- nrow(Tr)

    # Initialize the covariance matrix of T.
    Omega33 <- matrix(0, m3, m3)

    # If there are no triples, return the empty matrix.
    if (m3 == 0) return(Omega33)

    # Each entry of T is standardized, so the diagonal is 1.
    diag(Omega33) <- 1

    # Fill the off-diagonal covariances analytically.
    for (row in seq_len(m3)) {

      # Indices of the first three-way interaction.
      i <- Tr[row, 1]
      j <- Tr[row, 2]
      l <- Tr[row, 3]

      for (col in seq_len(m3)) {

        # The diagonal has already been set to 1.
        if (row == col) next

        # Indices of the second three-way interaction.
        a <- Tr[col, 1]
        b <- Tr[col, 2]
        c <- Tr[col, 3]

        # Fill Cov(T_ijl, T_abc).
        Omega33[row, col] <- cov_T3(i, j, l, a, b, c, Omega11)
      }
    }

    Omega33
  }

  # Build Omega13, the covariance matrix between the main-effect block C and
  # the standardized three-way interaction block T.
  build_Omega13 <- function(Omega11) {

    # List all triples that define entries of T.
    Tr <- all_triples(nrow(Omega11))

    # Number of three-way interaction terms.
    m3 <- nrow(Tr)

    # Initialize the k x m3 covariance matrix Cov(C, T).
    Omega13 <- matrix(0, k, m3)

    # If there are no triples, return the empty matrix.
    if (m3 == 0) return(Omega13)

    # Fill each entry Cov(C_a, T_ijl).
    for (col in seq_len(m3)) {

      # Indices of the three-way interaction.
      i <- Tr[col, 1]
      j <- Tr[col, 2]
      l <- Tr[col, 3]

      for (a in seq_len(k)) {
        Omega13[a, col] <- cov_C_T(a, i, j, l, Omega11)
      }
    }

    Omega13
  }

  # ---- solve scales ------------------------------------------------------------------------------
  #
  # This function finds the two scale factors used in Appendix A:
  #
  # - sL for the main-effect coefficient direction b_L;
  # - s  for the interaction coefficient directions b_2 and b_3.
  #
  # The unscaled quantities are:
  #
  # A = b_L' Omega11 b_L,      the main-effect variance;
  # B = b_2' Omega22 b_2,      the two-way interaction variance;
  # G = b_3' Omega33 b_3,      the three-way interaction variance;
  # D = b_L' Omega13 b_3,      the C/T covariance term.
  #
  # When D is non-zero, the cross term is split symmetrically between Vlin and
  # VNL, exactly as in Appendix A.

  solve_scales <- function(A, B, G, D, Vlin_star, VNL_star) {

    # If the target main-effect contribution is zero, the main-effect scale is
    # zero. Then only the interaction block is scaled.
    if (Vlin_star == 0) {

      # At least one interaction direction must have positive variance.
      if ((B + G) <= 0)
        stop("B + G must be > 0 when the interaction contribution is positive.")

      # With sL = 0, Var(eta_t) from interactions is s^2 (B + G).
      return(list(sL = 0, s = sqrt(VNL_star / (B + G)), r = Inf))
    }

    # If the target interaction contribution is zero, the interaction scale is
    # zero. Then only the main-effect block is scaled.
    if (VNL_star == 0) {

      # The sampled main-effect direction must have positive variance.
      if (A <= 0)
        stop("A must be > 0.")

      # With s = 0, Var(eta_t) from main effects is sL^2 A.
      return(list(sL = sqrt(Vlin_star / A), s = 0, r = 0))
    }

    # If no interaction block has variance and D is also zero, this reduces to
    # the linear case.
    if ((B + G) == 0 && D == 0) {

      # The sampled main-effect direction must have positive variance.
      if (A <= 0)
        stop("A must be > 0.")

      return(list(sL = sqrt(Vlin_star / A), s = 0, r = 0))
    }

    # If the main-effect and three-way blocks are uncorrelated, the split between
    # Vlin and VNL is direct.
    if (D == 0) {

      # Both relevant unscaled variance components must be positive.
      if (A <= 0)
        stop("A must be > 0.")
      if ((B + G) <= 0)
        stop("B + G must be > 0.")

      return(list(
        sL = sqrt(Vlin_star / A),
        s = sqrt(VNL_star / (B + G)),
        r = sqrt(VNL_star / (B + G)) / sqrt(Vlin_star / A)
      ))
    }

    # Otherwise, solve the quadratic equation in r = s / sL:
    #
    # Vlin_star (B + G) r^2 + D (Vlin_star - VNL_star) r - VNL_star A = 0.
    qa <- Vlin_star * (B + G)
    qb <- D * (Vlin_star - VNL_star)
    qc <- -VNL_star * A

    # Compute the discriminant. It should be non-negative under the assumptions
    # in Appendix A, apart from small numerical error.
    disc <- qb^2 - 4 * qa * qc

    # Reject genuinely negative discriminants.
    if (disc < 0)
      stop("No real solution for the ratio r = s / sL.")

    # Avoid numerical problems if disc is a tiny negative number rounded below 0.
    disc <- max(disc, 0)

    # Compute both roots of the quadratic.
    r1 <- (-qb + sqrt(disc)) / (2 * qa)
    r2 <- (-qb - sqrt(disc)) / (2 * qa)

    # A root is admissible if r is positive and gives A + rD > 0, because
    # sL^2 = Vlin_star / (A + rD).
    admissible <- function(r) is.finite(r) && (r > 0) && (A + r * D > 0)

    # Appendix A implies that there is exactly one positive admissible root when
    # Vlin_star > 0, VNL_star > 0, and B + G > 0.
    r <- if (admissible(r1)) r1 else if (admissible(r2)) r2 else NA_real_

    # Stop if numerical problems or incompatible inputs prevent a valid root.
    if (!is.finite(r))
      stop("No admissible positive solution for r = s / sL.")

    # Once r is known, compute sL and s.
    sL <- sqrt(Vlin_star / (A + r * D))
    s <- r * sL

    list(sL = sL, s = s, r = r)
  }

  # ---- covariance blocks -------------------------------------------------------------------------
  #
  # The covariance blocks depend only on Omega11 and on which interaction terms
  # are included. They are computed once and then reused during sampling.

  # Initialize the covariance blocks for Z, T, and Cov(C, T).
  Omega22 <- matrix(0, m2, m2)
  Omega33 <- matrix(0, m3, m3)
  Omega13 <- matrix(0, k, m3)

  # Compute the two-way interaction covariance block if Z is included.
  if (include_2way) Omega22 <- build_Omega22(Omega11)

  # Compute the three-way interaction covariance block and the C/T covariance
  # block if T is included.
  if (include_3way) {
    Omega33 <- build_Omega33(Omega11)
    Omega13 <- build_Omega13(Omega11)
  }

  # ---- feature names -----------------------------------------------------------------------------
  #
  # These names label the entries of M and the columns of the returned Delta_t.

  # Names for the main-effect block C.
  C_names <- paste0("C", seq_len(k))

  # Names for the two-way interaction block Z.
  Z_names <- if (include_2way) {
    combn(seq_len(k), 2, FUN = function(ix) paste0("Z", paste(ix, collapse = ":")))
  } else character(0)

  # Names for the three-way interaction block T.
  T_names <- if (include_3way) {
    combn(seq_len(k), 3, FUN = function(ix) paste0("T", paste(ix, collapse = ":")))
  } else character(0)

  # Full list of feature names in the same order as M = (C', Z', T')'.
  feature_names <- c(C_names, Z_names, T_names)

  # Attach names to Omega22 if Z exists.
  if (m2 > 0) {
    rownames(Omega22) <- Z_names
    colnames(Omega22) <- Z_names
  }

  # Attach names to Omega33 and Omega13 if T exists.
  if (m3 > 0) {
    rownames(Omega33) <- T_names
    colnames(Omega33) <- T_names
    rownames(Omega13) <- C_names
    colnames(Omega13) <- T_names
  }

  # ---- full Omega matrix -------------------------------------------------------------------------
  #
  # This helper builds the full covariance matrix Omega of M:
  #
  #           [ Omega11    0      Omega13 ]
  #   Omega = [    0    Omega22      0    ]
  #           [ Omega13'   0      Omega33 ]
  #
  # The zero blocks follow from the vanishing odd-order moments of centered
  # jointly Gaussian base confounders.

  build_full_Omega <- function(Omega11, Omega22, Omega33, Omega13) {

    # Total number of entries in M.
    p_total <- k + nrow(Omega22) + nrow(Omega33)

    # Initialize the full covariance matrix.
    Omega <- matrix(0, p_total, p_total)

    # Indices of the C, Z, and T blocks inside M.
    idx_C <- seq_len(k)
    idx_Z <- if (m2 > 0) (k + 1):(k + m2) else integer(0)
    idx_T <- if (m3 > 0) (k + m2 + 1):(k + m2 + m3) else integer(0)

    # Place Omega11 in the C/C block.
    Omega[idx_C, idx_C] <- Omega11

    # Place Omega22 in the Z/Z block.
    if (m2 > 0) Omega[idx_Z, idx_Z] <- Omega22

    # Place Omega33 in the T/T block.
    if (m3 > 0) Omega[idx_T, idx_T] <- Omega33

    # Place Omega13 and its transpose in the C/T and T/C blocks.
    if (m3 > 0) {
      Omega[idx_C, idx_T] <- Omega13
      Omega[idx_T, idx_C] <- t(Omega13)
    }

    # Attach feature names to the rows and columns.
    rownames(Omega) <- feature_names
    colnames(Omega) <- feature_names

    Omega
  }

  # ---- sample one row ----------------------------------------------------------------------------
  #
  # This helper samples one coefficient vector delta_t for one generic outcome.
  # The main function calls it twice: once for X and once for Y.

  sample_one_row <- function() {

    # Sample the initial coefficient direction for the main-effect block.
    b_L <- rnorm(k)

    # Sample the initial coefficient direction for the interaction block.
    b_int <- if (m_int > 0) rnorm(m_int) else numeric(0)

    # If requested, make all sampled directions non-negative before scaling.
    if (force_positive) {
      b_L <- abs(b_L)
      if (m_int > 0) b_int <- abs(b_int)
    }

    # Split the interaction direction into b_2 and b_3.
    b_2 <- if (m2 > 0) b_int[seq_len(m2)] else numeric(0)
    b_3 <- if (m3 > 0) b_int[(m2 + 1):(m2 + m3)] else numeric(0)

    # A is the unscaled main-effect variance b_L' Omega11 b_L.
    A <- as.numeric(t(b_L) %*% Omega11 %*% b_L)

    # B is the unscaled two-way interaction variance b_2' Omega22 b_2.
    B <- if (m2 > 0) {
      as.numeric(t(b_2) %*% Omega22 %*% b_2)
    } else 0

    # G is the unscaled three-way interaction variance b_3' Omega33 b_3.
    G <- if (m3 > 0) {
      as.numeric(t(b_3) %*% Omega33 %*% b_3)
    } else 0

    # D is the unscaled covariance term b_L' Omega13 b_3.
    D <- if (m3 > 0) {
      as.numeric(t(b_L) %*% Omega13 %*% b_3)
    } else 0

    # If M contains only C, this is the linear case from Appendix A. A single
    # scale factor sL is sufficient.
    if (m_int == 0) {
      if (A <= 0) stop("A must be > 0; check Omega11.")
      sL <- sqrt(Vlin_star / A)
      s <- 0
    } else {

      # Otherwise, solve for sL and s using A, B, G, and D.
      scales <- solve_scales(
        A = A,
        B = B,
        G = G,
        D = D,
        Vlin_star = Vlin_star,
        VNL_star = VNL_star
      )

      sL <- scales$sL
      s <- scales$s
    }

    # Apply the scale factors to obtain the final coefficient blocks.
    delta_L <- sL * b_L
    delta_2 <- if (m2 > 0) s * b_2 else numeric(0)
    delta_3 <- if (m3 > 0) s * b_3 else numeric(0)

    # Return delta_t in the same order as M = (C', Z', T')'.
    c(delta_L, delta_2, delta_3)
  }

  # ---- main sampling loop ------------------------------------------------------------------------
  #
  # The procedure samples one row for X and one row for Y. If coefficient bounds
  # are imposed through min_abs and max_abs, sampling is repeated until both rows
  # satisfy those bounds or max_tries is reached.

  for (attempt in seq_len(max_tries)) {

    # Sample the coefficient vector for the direct contribution to X_t.
    delta_X <- sample_one_row()

    # Sample the coefficient vector for the direct contribution to Y_t.
    delta_Y <- sample_one_row()

    # Check whether all sampled coefficients satisfy the requested bounds.
    if (all(abs(c(delta_X, delta_Y)) >= min_abs &
            abs(c(delta_X, delta_Y)) <= max_abs)) {

      # Stack the two coefficient vectors into Delta_t.
      Delta_t <- rbind(delta_X, delta_Y)
      rownames(Delta_t) <- c("X", "Y")
      colnames(Delta_t) <- feature_names

      # Build the full covariance matrix Omega of M.
      Omega <- build_full_Omega(Omega11, Omega22, Omega33, Omega13)

      # Return the sampled coefficient matrix and the full feature covariance
      # matrix. These are sufficient to verify Var(eta_t) = R2_total.
      return(list(
        Delta_t = Delta_t,
        Delta = Delta_t,
        Omega = Omega,
        Omega11 = Omega11,
        Omega22 = Omega22,
        Omega33 = Omega33,
        Omega13 = Omega13
      ))
    }
  }

  stop("Failed to sample Delta_t within max_tries. Try relaxing bounds or increasing max_tries.")
}

