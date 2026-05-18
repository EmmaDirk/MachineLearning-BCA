# =================================================================================================
#
# This function simulates panel data under a cross-lagged data-generating process
# with time-invariant confounders whose effects may vary over time.
#
# The goal is to generate two observed variables, X and Y, measured over T waves
# for N individuals, where X and Y are confounded by a time-invariant set of
# confounder features C.
#
# The simulation enforces that the observed covariance matrix of (X_t, Y_t)
# equals Sigma at every wave:
#
#   Var(X_t) = 1
#   Var(Y_t) = 1
#   Cov(X_t, Y_t) = Sigma[1,2]
#
# Naming logic:
#
# Phi        = lag matrix of autoregressive and cross-lagged effects
# Sigma      = target covariance matrix of (X_t, Y_t)
# Omega11    = covariance matrix of base confounders C1..Ck
# Delta_list = list of confounder effect matrices Delta_t
#
# Delta_t maps confounder features into (X_t, Y_t):
#
#   W_t = Phi W_{t-1} + Delta_t C + epsilon_t
#
# where
#
#   W_t       = (X_t, Y_t)'
#   C         = vector of confounder features
#   epsilon_t ~ (0, Psi_t)
#
# Because the same confounders affect every wave, the lagged term Phi W_{t-1}
# and the direct confounder term Delta_t C are correlated. Therefore we track:
#
#   M_t = Cov(W_t, C)
#
# with recursion:
#
#   M_t = Phi M_{t-1} + Delta_t Omega
#
# The innovation covariance is chosen such that Var(W_t) = Sigma:
#
#   Psi_t =
#     Sigma
#     - Phi Sigma Phi'
#     - Delta_t Omega Delta_t'
#     - Phi M_{t-1} Delta_t'
#     - Delta_t M_{t-1}' Phi'
#
# This guarantees the requested observed covariance structure at every wave.
# =================================================================================================


# ---- true confounder R2 trajectory ---------------------------------------------------------------

# This helper computes the true confounder R^2 trajectory analytically.
#
# Why this helper exists:
# - the simulation study already stores out-of-fold R^2 from the residualisers
# - however, those are model-based prediction diagnostics
# - for interpretation, we also want the population benchmark:
#   the true fraction of variance in X_t and Y_t that is attributable to the
#   confounder feature vector C at each wave
#
# In this DGM, the total confounder effect at wave t is not just the direct
# contemporaneous term Delta_t C. It also includes all earlier confounder effects
# that have propagated through the lagged system Phi.
#
# Writing
#
#   W_t = Phi W_{t-1} + Delta_t C + epsilon_t
#
# the total linear mapping from C into W_t can be written as
#
#   B_1 = Delta_1
#   B_t = Phi B_{t-1} + Delta_t
#
# so that
#
#   E(W_t | C) = B_t C .
#
# Therefore the population covariance explained by the confounders at wave t is
#
#   Var(E(W_t | C)) = B_t Omega B_t'
#
# where Omega is the full population covariance matrix of the confounder feature
# vector C, including interaction features in the same order as the Delta columns.
#
# Because this simulation enforces Var(X_t) = Var(Y_t) = 1 through Sigma, the
# diagonal elements of B_t Omega B_t' are already the true wave-specific R^2
# values for X_t and Y_t.

compute_true_confounder_r2 <- function(
    T,                                                            # number of observed waves to keep
    Phi,                                                          # lag matrix
    Delta_list,                                                   # full internal Delta trajectory, including burn-in
    Omega11,                                                      # covariance matrix of base confounders
    Sigma,                                                        # target covariance matrix of (X_t, Y_t)
    burn_in = 0L,                                                 # number of internal burn-in waves
    eig_tol = 1e-10                                               # tolerance for positive semidefiniteness
){

  # ---- input checks ------------------------------------------------------------------------------

  # We repeat the key checks here because this helper is intended to be callable
  # directly from the replication wrapper, without assuming that simulate_panel_data()
  # has already been called.

  # Check that T is a positive integer.

  if (!is.numeric(T) || length(T) != 1 || is.na(T) || T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # Check that burn_in is a non-negative integer.

  if (!is.numeric(burn_in) || length(burn_in) != 1 || is.na(burn_in) ||
      burn_in < 0 || burn_in != as.integer(burn_in))
    stop("burn_in must be a non-negative integer.")

  # Check that Phi and Sigma are 2x2 matrices.

  if (!is.matrix(Phi) || !all(dim(Phi) == c(2,2)))
    stop("Phi must be a 2x2 matrix.")
  if (!is.matrix(Sigma) || !all(dim(Sigma) == c(2,2)))
    stop("Sigma must be a 2x2 matrix.")

  # Sigma is the marginal variance target for (X_t, Y_t).

  if (!isTRUE(all.equal(Sigma, t(Sigma))))
    stop("Sigma must be symmetric.")

  # In this study the diagonal is fixed to 1, so the diagonal of the explained
  # covariance matrix is directly interpretable as R^2.

  if (!isTRUE(all.equal(diag(Sigma), c(1, 1))))
    stop("Sigma must have 1 on the diagonal.")

  # Omega11 must describe standardized Gaussian base confounders.

  if (!is.matrix(Omega11))
    stop("Omega11 must be a matrix.")
  if (!isTRUE(all.equal(Omega11, t(Omega11))))
    stop("Omega11 must be symmetric.")
  if (!isTRUE(all.equal(diag(Omega11), rep(1, nrow(Omega11)))))
    stop("Base confounders must be standardized (diag(Omega11)=1).")

  # Delta_list must contain the full internal trajectory, including burn-in.

  if (!is.list(Delta_list) || length(Delta_list) == 0)
    stop("Delta_list must be a non-empty list.")

  # Coerce counts once.

  T_obs <- as.integer(T)
  burn_in <- as.integer(burn_in)
  T_total <- T_obs + burn_in

  # The trajectory generators in 02_delta_trajectory.R already produce the full
  # internal path, so we require that exact convention here too.

  if (length(Delta_list) != T_total) {
    stop(
      "Delta_list must have length T + burn_in. ",
      "You supplied length ", length(Delta_list),
      ", but expected ", T_total, "."
    )
  }

  # All Delta matrices must share the same feature names and shape.

  feature_names <- colnames(Delta_list[[1]])
  if (is.null(feature_names))
    stop("Delta matrices must have column names.")

  # Check that every Delta matrix has the same dimensions and same feature order.

  for (t in seq_along(Delta_list)) {

    # Every element must be a matrix.

    if (!is.matrix(Delta_list[[t]]))
      stop("Every entry in Delta_list must be a matrix.")

    # All Delta matrices must be 2 x p.

    if (!all(dim(Delta_list[[t]]) == dim(Delta_list[[1]])))
      stop("All Delta matrices in Delta_list must have identical dimensions.")

    # The feature order must remain the same over time.

    if (!identical(colnames(Delta_list[[t]]), feature_names))
      stop("All Delta matrices in Delta_list must have the same column names in the same order.")
  }

  # ---- positive semidefiniteness check -----------------------------------------------------------

  # This is used as a safety check for derived covariance matrices.

  check_psd <- function(S) {

    # Enforce symmetry numerically before checking the eigenvalues.

    S <- (S + t(S)) / 2

    # Eigen decomposition.

    evd <- eigen(S, symmetric = TRUE)

    # Extract eigenvalues.

    vals <- evd$values

    # Stop if any eigenvalue is meaningfully negative.

    if (any(vals < -eig_tol))
      stop("A covariance matrix is not positive semidefinite.")

    # Otherwise return the symmetrized matrix.

    S
  }

  # ---- feature-name parsing ----------------------------------------------------------------------

  # This helper translates names such as:
  # - C1       -> 1
  # - Z1:2     -> 1,2
  # - T1:2:3   -> 1,2,3
  #
  # This is the naming convention used in the simulation code.

  parse_feature <- function(name) {

    # Remove the Appendix A block prefix from the first component and keep the numeric indices.

    parts <- strsplit(name, ":", fixed = TRUE)[[1]]
    parts[1] <- sub("^[CcZzTt]", "", parts[1])

    as.integer(parts)
  }

  # ---- Delta vector storage ----------------------------------------------------------------------

  # This helper flattens one Delta_t into a named vector.
  #
  # The user asked to keep the true Delta information in the saved output, but not
  # necessarily as one separate scalar column per coefficient.
  #
  # To support that workflow, we flatten each 2 x p Delta_t matrix into one named
  # numeric vector. This vector is then easy to store in a single list-column of the
  # final results data frame, while still preserving the full coefficient information.
  #
  # Naming convention:
  # - coefficients for the X equation are prefixed with "x__"
  # - coefficients for the Y equation are prefixed with "y__"
  # - the original confounder-feature names are retained after the prefix
  #
  # Example:
  #
  #   Delta_t[, c("C1", "Z1:2")]
  #
  # becomes a named vector with entries such as
  #
  #   x__C1, x__Z1:2, y__C1, y__Z1:2
  #
  # We flatten row-wise so that the X-row coefficients appear first, followed by the
  # Y-row coefficients.

  flatten_Delta_matrix <- function(Delta_t) {

    # Require a 2 x p matrix because the DGM has exactly two observed variables.

    if (!is.matrix(Delta_t) || nrow(Delta_t) != 2L)
      stop("Delta_t must be a 2 x p matrix.")

    # Require feature names so the stored vector remains interpretable later.

    if (is.null(colnames(Delta_t)))
      stop("Delta_t must have column names.")

    # Build readable names for both outcome rows.

    coef_names <- c(
      paste0("x__", colnames(Delta_t)),
      paste0("y__", colnames(Delta_t))
    )

    # Flatten row-wise: first all X coefficients, then all Y coefficients.

    out <- c(Delta_t[1, ], Delta_t[2, ])
    names(out) <- coef_names

    # Return one named numeric vector.

    out
  }

  # ---- analytic covariance of confounder features ------------------------------------------------

  # We intentionally compute the feature covariance matrix analytically instead of
  # estimating cov(C_full) from a finite sample. This gives the population Omega
  # that the DGM is based on, and therefore the population benchmark for the true
  # confounder-explained variance.

  # For centered Gaussian variables, the sixth moment equals the sum over all 15
  # pairings of the 6 indices. We build those pairings once and reuse them.

  pairings6 <- local({

    rec <- function(v) {

      # If there are no indices left, return one empty pairing structure.

      if (length(v) == 0) return(list(list()))

      # If there are exactly two indices left, they form the final pair.

      if (length(v) == 2) return(list(list(c(v[1], v[2]))))

      # Otherwise recursively pair the first element with each later element.

      first <- v[1]
      out <- list()

      for (m in 2:length(v)) {

        # Remove the chosen pair and recurse on the remainder.

        rest <- v[-c(1, m)]
        sub  <- rec(rest)

        # Append the chosen pair to each recursive solution.

        for (s in sub) {
          out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
        }
      }

      out
    }

    rec(1:6)
  })

  # Compute E[X1 X2 X3 X4 X5 X6] for centered Gaussian variables using
  # Isserlis' theorem: sum over all pairings of products of covariances.

  sixth_moment_gaussian <- function(idx, Omega11) {

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

  # Variance of the raw two-way product C_i * C_l.
  # For standardized Gaussian variables this is 1 + rho^2.

  var_raw_2way <- function(i, l, Omega11) {

    1 + Omega11[i,l]^2
  }

  # Variance of the raw three-way product C_i * C_l * C_m.
  # This is the denominator used to standardize the three-way interaction.

  var_raw_3way <- function(i, l, m, Omega11) {

    rho_il <- Omega11[i,l]
    rho_im <- Omega11[i,m]
    rho_lm <- Omega11[l,m]

    1 +
      2*rho_il^2 +
      2*rho_im^2 +
      2*rho_lm^2 +
      8*rho_il*rho_im*rho_lm
  }

  # Covariance between two main effects.

  cov_main_main <- function(a, b, Omega11) {

    Omega11[a,b]
  }

  # Covariance between a main effect and a standardized two-way interaction.
  # For centered Gaussian variables this is 0.

  cov_main_2way <- function(a, i, l, Omega11) {

    0
  }

  # Covariance between a main effect and a standardized three-way interaction.

  cov_main_3way <- function(a, i, l, m, Omega11) {

    num <- Omega11[a,i] * Omega11[l,m] +
           Omega11[a,l] * Omega11[i,m] +
           Omega11[a,m] * Omega11[i,l]

    den <- sqrt(var_raw_3way(i, l, m, Omega11))

    num / den
  }

  # Covariance between two standardized two-way interactions.

  cov_2way_2way <- function(i, l, p, q, Omega11) {

    num <- Omega11[i,p] * Omega11[l,q] +
           Omega11[i,q] * Omega11[l,p]

    den <- sqrt(
      var_raw_2way(i, l, Omega11) *
      var_raw_2way(p, q, Omega11)
    )

    num / den
  }

  # Covariance between a standardized two-way interaction and a standardized
  # three-way interaction. For centered Gaussian variables this is 0.

  cov_2way_3way <- function(i, l, p, q, r, Omega11) {

    0
  }

  # Covariance between two standardized three-way interactions.

  cov_3way_3way <- function(i, l, m, p, q, r, Omega11) {

    num <- sixth_moment_gaussian(c(i, l, m, p, q, r), Omega11)

    den <- sqrt(
      var_raw_3way(i, l, m, Omega11) *
      var_raw_3way(p, q, r, Omega11)
    )

    num / den
  }

  # Determine what kind of feature we are dealing with.

  feature_type <- function(ix) {

    order <- length(ix)

    if (order == 1) return("main")
    if (order == 2) return("2way")
    if (order == 3) return("3way")

    stop("Only up to three-way interactions supported.")
  }

  # Compute the covariance between any two supported features.

  cov_feature_pair <- function(ix1, ix2, Omega11) {

    type1 <- feature_type(ix1)
    type2 <- feature_type(ix2)

    # Main with main.

    if (type1 == "main" && type2 == "main") {
      return(cov_main_main(ix1[1], ix2[1], Omega11))
    }

    # Main with two-way.

    if (type1 == "main" && type2 == "2way") {
      return(cov_main_2way(ix1[1], ix2[1], ix2[2], Omega11))
    }

    # Two-way with main.

    if (type1 == "2way" && type2 == "main") {
      return(cov_main_2way(ix2[1], ix1[1], ix1[2], Omega11))
    }

    # Main with three-way.

    if (type1 == "main" && type2 == "3way") {
      return(cov_main_3way(ix1[1], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # Three-way with main.

    if (type1 == "3way" && type2 == "main") {
      return(cov_main_3way(ix2[1], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # Two-way with two-way.

    if (type1 == "2way" && type2 == "2way") {
      return(cov_2way_2way(ix1[1], ix1[2], ix2[1], ix2[2], Omega11))
    }

    # Two-way with three-way.

    if (type1 == "2way" && type2 == "3way") {
      return(cov_2way_3way(ix1[1], ix1[2], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # Three-way with two-way.

    if (type1 == "3way" && type2 == "2way") {
      return(cov_2way_3way(ix2[1], ix2[2], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # Three-way with three-way.

    if (type1 == "3way" && type2 == "3way") {
      return(cov_3way_3way(ix1[1], ix1[2], ix1[3], ix2[1], ix2[2], ix2[3], Omega11))
    }

    stop("Unsupported feature combination.")
  }

  # Build the full covariance matrix of the confounder feature vector in exactly
  # the same order as feature_names.

  build_full_Omega_from_features <- function(feature_names, Omega11) {

    # p is the total number of requested features.

    p <- length(feature_names)

    # Parse every feature name once.

    parsed_features <- lapply(feature_names, parse_feature)

    # Initialize Omega.

    Omega <- matrix(0, p, p)
    rownames(Omega) <- feature_names
    colnames(Omega) <- feature_names

    # Fill the upper triangle and mirror it to the lower triangle.

    for (a in seq_len(p)) {
      for (b in a:p) {

        val <- cov_feature_pair(parsed_features[[a]], parsed_features[[b]], Omega11)

        Omega[a,b] <- val
        Omega[b,a] <- val
      }
    }

    # Numerically symmetrize and check PSD.

    check_psd(Omega)
  }

  # ---- build population feature covariance -------------------------------------------------------

  # Omega_full is the population covariance matrix of the confounder feature vector
  # in the exact same column order used by every Delta_t.

  Omega_full <- build_full_Omega_from_features(feature_names, Omega11)

  # ---- propagate total confounder effects through time -------------------------------------------

  # B_t collects the total mapping from C into W_t after both:
  # - the direct same-wave confounder effect Delta_t
  # - all earlier confounder effects that were carried forward by Phi
  #
  # Recursion:
  #   B_1 = Delta_1
  #   B_t = Phi B_{t-1} + Delta_t

  B_list <- vector("list", T_total)

  # Baseline wave.

  B_list[[1]] <- Delta_list[[1]]

  # Later waves.

  if (T_total >= 2L) {
    for (t in 2:T_total) {
      B_list[[t]] <- Phi %*% B_list[[t - 1]] + Delta_list[[t]]
    }
  }

  # ---- compute true explained covariance by wave -------------------------------------------------

  # For each internal wave t:
  #
  #   Var(E(W_t | C)) = B_t Omega_full B_t'
  #
  # The diagonal entries are the population confounder-explained variances in X_t
  # and Y_t. Because Sigma has unit diagonal, those are also the true R^2 values.

  out <- data.frame(
    t_internal = seq_len(T_total),
    wave_observed = rep(NA_integer_, T_total),
    true_r2_x = rep(NA_real_, T_total),
    true_r2_y = rep(NA_real_, T_total),
    true_cov_xy_confounder = rep(NA_real_, T_total),
    stringsAsFactors = FALSE
  )

  # Prepare one list-column that stores the direct true Delta_t coefficients for
  # every internal wave. Each entry will be a named numeric vector created by the
  # helper above.

  out$true_delta_t_vector <- vector("list", T_total)

  # Internal indices that survive burn-in and are therefore analyzed.

  keep_idx <- seq.int(from = burn_in + 1L, to = T_total)

  # Label the analyzed waves on the observed 1..T scale.

  out$wave_observed[keep_idx] <- seq_len(T_obs)

  # Compute the confounder-explained covariance wave by wave.

  for (t in seq_len(T_total)) {

    # Total covariance explained by the confounders at wave t.

    V_conf_t <- B_list[[t]] %*% Omega_full %*% t(B_list[[t]])

    # Enforce symmetry / PSD numerically.

    V_conf_t <- check_psd(V_conf_t)

    # Because Var(X_t)=Var(Y_t)=1 by construction, the diagonal entries already
    # equal the true R^2 values. We still divide by Sigma's diagonal explicitly
    # to keep the formula transparent and robust to future extensions.

    out$true_r2_x[t] <- V_conf_t[1,1] / Sigma[1,1]
    out$true_r2_y[t] <- V_conf_t[2,2] / Sigma[2,2]

    # Store the confounder-explained covariance between X_t and Y_t as well.
    # This was not requested as a main output, but it is often useful for diagnostics.

    out$true_cov_xy_confounder[t] <- V_conf_t[1,2]

    # Also store the direct Delta_t coefficients themselves for this wave.
    # This makes it possible to inspect the exact confounder-effect matrix later
    # from the saved results data frame, without creating one scalar column per
    # coefficient.

    out$true_delta_t_vector[[t]] <- flatten_Delta_matrix(Delta_list[[t]])
  }

  # Mark the Delta vector column explicitly as a list-column. This prevents
  # data.frame from trying to simplify the nested vectors when the object is later
  # merged into the final simulation output.

  out$true_delta_t_vector <- I(out$true_delta_t_vector)

  # Keep only the analyzed waves in a separate compact object as well.

  out_observed <- out[keep_idx, , drop = FALSE]

  # Relabel rows cleanly.

  rownames(out) <- NULL
  rownames(out_observed) <- NULL

  # Return both versions:
  # - all_waves: the full internal trajectory, including burn-in
  # - observed_waves: only the analyzed waves 1..T
  #
  # Both data frames now also include a list-column named true_delta_t_vector that
  # stores the direct Delta_t matrix for each wave in flattened named-vector form.

  list(
    all_waves = out,
    observed_waves = out_observed,
    Omega_full = Omega_full,
    B_list = B_list
  )
}


# ---- panel-data simulation -----------------------------------------------------------------------

simulate_panel_data <- function(
    N,                                                            # number of individuals
    T,                                                            # number of observed waves to keep
    Phi,                                                          # lag matrix (diagonal elements are autoregressive effects)
    Delta_list,                                                   # list of Delta matrices
    Omega11,                                                      # covariance matrix of base confounders
    Sigma,                                                        # desired covariance matrix of (X_t, Y_t)
    burn_in = 0L,                                                 # number of initial waves to discard
    seed = NULL,                                                  # seed for reproducibility
    eig_tol = 1e-10                                               # tolerance for positive semidefiniteness
){

  # ---- input checks ------------------------------------------------------------------------------

  # Check that N and T are positive integers.

  if (!is.numeric(N) || N < 1 || N != as.integer(N))
    stop("N must be a positive integer.")
  if (!is.numeric(T) || T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")
  if (!is.numeric(burn_in) || burn_in < 0 || burn_in != as.integer(burn_in))
    stop("burn_in must be a non-negative integer.")

  # Observed and total number of waves.

  T_obs <- as.integer(T)
  burn_in <- as.integer(burn_in)
  T_total <- T_obs + burn_in

  # Check that Phi and Sigma are 2x2 matrices.

  if (!is.matrix(Phi) || !all(dim(Phi) == c(2,2)))
    stop("Phi must be a 2x2 matrix.")
  if (!is.matrix(Sigma) || !all(dim(Sigma) == c(2,2)))
    stop("Sigma must be a 2x2 matrix.")

  # Check that Sigma is symmetric.

  if (!isTRUE(all.equal(Sigma, t(Sigma))))
    stop("Sigma must be symmetric.")
  
  # Check that Sigma has 1 on the diagonal.

  if (!isTRUE(all.equal(diag(Sigma), c(1,1))))
    stop("Sigma must have 1 on the diagonal.")

  # Check that Omega11 is a square matrix.

  if (!is.matrix(Omega11))
    stop("Omega11 must be a matrix.")

  # Check that Omega11 is symmetric.

  if (!isTRUE(all.equal(Omega11, t(Omega11))))
    stop("Omega11 must be symmetric.")

  # Check that Omega11 has 1 on the diagonal.

  if (!isTRUE(all.equal(diag(Omega11), rep(1, nrow(Omega11)))))
    stop("Base confounders must be standardized (diag(Omega11)=1).")

  # Check that Delta_list is a list.

  if (!is.list(Delta_list))
    stop("Delta_list must be a list.")

  # Delta_list must follow one strict convention:
  # it must represent the full internal trajectory, including burn-in waves.
  # The trajectory generators in 02_delta_trajectory.R already handle this.

  if (length(Delta_list) != T_total) {
    stop(
      "Delta_list must have length T + burn_in. ",
      "You supplied length ", length(Delta_list),
      ", but expected ", T_total, "."
    )
  }

  # Check that Delta matrices have column names.

  feature_names <- colnames(Delta_list[[1]])
  if (is.null(feature_names))
    stop("Delta matrices must have column names.")

  # Keep any names supplied by the trajectory generators. If Delta_list has no
  # names, fall back to positional names.

  if (is.null(names(Delta_list))) {
    names(Delta_list) <- paste0("t", seq_len(T_total))
  }

  # k is the number of base confounders.

  k <- nrow(Omega11)

  # Use Appendix A main-effect names of the form C1, ..., Ck.

  lin_names <- paste0("C", 1:k)

  # ---- positive semidefiniteness check -----------------------------------------------------------

  # Check that covariance matrices are positive semidefinite.

  check_psd <- function(S) {

    # Symmetrize the matrix numerically.

    S <- (S + t(S))/2

    # Eigen decomposition.

    evd <- eigen(S, symmetric = TRUE)

    # Extract eigenvalues.

    vals <- evd$values

    # Stop if any eigenvalue is meaningfully negative.

    if (any(vals < -eig_tol))
      stop("A covariance matrix is not positive semidefinite.")

    # Otherwise return the symmetrized matrix.

    S
  }

  # ---- feature-name parsing ----------------------------------------------------------------------

  # This helper converts names such as C1, Z1:2, and T1:2:3 to their integer indices.

  parse_feature <- function(name) {

    # Remove the Appendix A block prefix from the first component and keep the numeric indices.

    parts <- strsplit(name, ":", fixed = TRUE)[[1]]
    parts[1] <- sub("^[CcZzTt]", "", parts[1])

    as.integer(parts)
  }

  # ---- analytic covariance of confounder features ------------------------------------------------

  # We do not need to estimate Cov(C_full) empirically from the simulated sample.
  # Instead, because all features are either:
  # - standardized main effects
  # - standardized two-way interactions
  # - standardized three-way interactions
  #
  # and because the base confounders are multivariate normal with covariance Omega11,
  # we can derive the full covariance matrix analytically.
  #
  # This gives the population covariance of the feature vector C, which is more
  # stable than cov(C_full) and does not depend on N.

  # For centered Gaussian variables, the sixth moment equals the sum over all 15
  # pairings of the 6 indices. We build those 15 pairings once and reuse them.

  pairings6 <- local({

    rec <- function(v) {

      # If there are no indices left, return one empty pairing structure.

      if (length(v) == 0) return(list(list()))

      # If there are exactly two indices left, they form the final pair.

      if (length(v) == 2) return(list(list(c(v[1], v[2]))))

      # Otherwise recursively pair the first element with each later element.

      first <- v[1]
      out <- list()

      for (m in 2:length(v)) {

        # Remove the chosen pair and recurse on the remainder.

        rest <- v[-c(1, m)]
        sub  <- rec(rest)

        # Append the chosen pair to each recursive solution.

        for (s in sub) {
          out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
        }
      }

      out
    }

    rec(1:6)
  })

  # Compute E[X1 X2 X3 X4 X5 X6] for centered Gaussian variables using
  # Isserlis' theorem: sum over all pairings of products of covariances.

  sixth_moment_gaussian <- function(idx, Omega11) {

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

  # Variance of the raw two-way product C_i * C_l.
  # For standardized Gaussian variables this is 1 + rho^2.

  var_raw_2way <- function(i, l, Omega11) {

    1 + Omega11[i,l]^2
  }

  # Variance of the raw three-way product C_i * C_l * C_m.
  # This is the denominator used to standardize the three-way interaction.

  var_raw_3way <- function(i, l, m, Omega11) {

    rho_il <- Omega11[i,l]
    rho_im <- Omega11[i,m]
    rho_lm <- Omega11[l,m]

    1 +
      2*rho_il^2 +
      2*rho_im^2 +
      2*rho_lm^2 +
      8*rho_il*rho_im*rho_lm
  }

  # Covariance between two main effects.

  cov_main_main <- function(a, b, Omega11) {

    Omega11[a,b]
  }

  # Covariance between a main effect and a standardized two-way interaction.
  # For centered Gaussian variables this is 0.

  cov_main_2way <- function(a, i, l, Omega11) {

    0
  }

  # Covariance between a main effect and a standardized three-way interaction.

  cov_main_3way <- function(a, i, l, m, Omega11) {

    num <- Omega11[a,i] * Omega11[l,m] +
           Omega11[a,l] * Omega11[i,m] +
           Omega11[a,m] * Omega11[i,l]

    den <- sqrt(var_raw_3way(i, l, m, Omega11))

    num / den
  }

  # Covariance between two standardized two-way interactions.

  cov_2way_2way <- function(i, l, p, q, Omega11) {

    num <- Omega11[i,p] * Omega11[l,q] +
           Omega11[i,q] * Omega11[l,p]

    den <- sqrt(
      var_raw_2way(i, l, Omega11) *
      var_raw_2way(p, q, Omega11)
    )

    num / den
  }

  # Covariance between a standardized two-way interaction and a standardized
  # three-way interaction. For centered Gaussian variables this is 0.

  cov_2way_3way <- function(i, l, p, q, r, Omega11) {

    0
  }

  # Covariance between two standardized three-way interactions.

  cov_3way_3way <- function(i, l, m, p, q, r, Omega11) {

    num <- sixth_moment_gaussian(c(i, l, m, p, q, r), Omega11)

    den <- sqrt(
      var_raw_3way(i, l, m, Omega11) *
      var_raw_3way(p, q, r, Omega11)
    )

    num / den
  }

  # Determine what kind of feature we are dealing with.

  feature_type <- function(ix) {

    order <- length(ix)

    if (order == 1) return("main")
    if (order == 2) return("2way")
    if (order == 3) return("3way")

    stop("Only up to three-way interactions supported.")
  }

  # Compute the covariance between any two supported features.

  cov_feature_pair <- function(ix1, ix2, Omega11) {

    type1 <- feature_type(ix1)
    type2 <- feature_type(ix2)

    # Main with main.

    if (type1 == "main" && type2 == "main") {
      return(cov_main_main(ix1[1], ix2[1], Omega11))
    }

    # Main with two-way.

    if (type1 == "main" && type2 == "2way") {
      return(cov_main_2way(ix1[1], ix2[1], ix2[2], Omega11))
    }

    # Two-way with main.

    if (type1 == "2way" && type2 == "main") {
      return(cov_main_2way(ix2[1], ix1[1], ix1[2], Omega11))
    }

    # Main with three-way.

    if (type1 == "main" && type2 == "3way") {
      return(cov_main_3way(ix1[1], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # Three-way with main.

    if (type1 == "3way" && type2 == "main") {
      return(cov_main_3way(ix2[1], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # Two-way with two-way.

    if (type1 == "2way" && type2 == "2way") {
      return(cov_2way_2way(ix1[1], ix1[2], ix2[1], ix2[2], Omega11))
    }

    # Two-way with three-way.

    if (type1 == "2way" && type2 == "3way") {
      return(cov_2way_3way(ix1[1], ix1[2], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # Three-way with two-way.

    if (type1 == "3way" && type2 == "2way") {
      return(cov_2way_3way(ix2[1], ix2[2], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # Three-way with three-way.

    if (type1 == "3way" && type2 == "3way") {
      return(cov_3way_3way(ix1[1], ix1[2], ix1[3], ix2[1], ix2[2], ix2[3], Omega11))
    }

    stop("Unsupported feature combination.")
  }

  # Build the full covariance matrix of the confounder feature vector in exactly
  # the same order as feature_names.

  build_full_Omega_from_features <- function(feature_names, Omega11) {

    # p is the total number of requested features.

    p <- length(feature_names)

    # Parse every feature name once.

    parsed_features <- lapply(feature_names, parse_feature)

    # Initialize Omega.

    Omega <- matrix(0, p, p)
    rownames(Omega) <- feature_names
    colnames(Omega) <- feature_names

    # Fill the upper triangle and mirror it to the lower triangle.

    for (a in seq_len(p)) {
      for (b in a:p) {

        val <- cov_feature_pair(parsed_features[[a]], parsed_features[[b]], Omega11)

        Omega[a,b] <- val
        Omega[b,a] <- val
      }
    }

    Omega
  }

  # ---- construct confounder feature matrix -------------------------------------------------------

  # Base confounders are simulated directly.
  # Interaction terms are constructed analytically and standardized using the same
  # formulas assumed in sample_delta_t().

  build_confounder_features <- function(C_base, feature_names, Omega11) {

    # N is the number of individuals.

    N <- nrow(C_base)

    # Initialize the matrix of confounder features.

    out <- matrix(NA, N, length(feature_names))
    colnames(out) <- feature_names

    # Fill every entry in the feature matrix.

    for (j in seq_along(feature_names)) {

      # Get the name of the feature.

      name <- feature_names[j]

      # Parse to get the interaction order.

      ix <- parse_feature(name)

      # The order is now the number of elements.

      order <- length(ix)

      # Main effect.

      if (order == 1) {

        # If order is 1, simply return the base confounder.

        out[,j] <- C_base[,ix]

      # Two-way interaction.
      } else if (order == 2) {

        # Make the index.

        i <- ix[1]; l <- ix[2]

        # Extract the interaction between the base confounders.
        # E[C_iC_l].

        rho <- Omega11[i,l]

        # Compute the raw interaction.

        raw <- C_base[,i] * C_base[,l]

        # Standardize using the formula.

        out[,j] <- (raw - rho) / sqrt(1 + rho^2)

      # Three-way interaction.
      } else if (order == 3) {

        # Make the index.

        i <- ix[1]; l <- ix[2]; m <- ix[3]

        # Get the covariances: E[C_iC_l], E[C_iC_m], E[C_lC_m].

        rho_il <- Omega11[i,l]
        rho_im <- Omega11[i,m]
        rho_lm <- Omega11[l,m]

        # Make the raw interaction.

        raw <- C_base[,i] * C_base[,l] * C_base[,m]

        # Compute the variance of that raw interaction.

        denom <- sqrt(
          1 +
          2*rho_il^2 +
          2*rho_im^2 +
          2*rho_lm^2 +
          8*rho_il*rho_im*rho_lm
        )

        # Standardize the output.

        out[,j] <- raw / denom

      } else {

        stop("Only up to three-way interactions supported.")
      }
    }

    out
  }

  # ---- seed --------------------------------------------------------------------------------------

  if (!is.null(seed))
    set.seed(seed)

  # ---- simulate base confounders -----------------------------------------------------------------

  C_base <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(0,k),
    sigma = Omega11
  )

  colnames(C_base) <- lin_names

  # ---- build full confounder feature matrix ------------------------------------------------------

  C_full <- build_confounder_features(
    C_base,
    feature_names,
    Omega11
  )

  # Analytic covariance of the confounder feature matrix.

  Omega_full <- build_full_Omega_from_features(
    feature_names,
    Omega11
  )

  # ---- containers --------------------------------------------------------------------------------

  Psi_list <- vector("list", T_total)
  M_list <- vector("list", T_total)
  W_list <- vector("list", T_total)

  # ---- wave 1 ------------------------------------------------------------------------------------

  # Extract the initial Delta_t matrix.

  Delta_initial <- Delta_list[[1]]

  # Psi1 = Sigma - Delta_initial %*% Omega %*% t(Delta_initial), since the
  # previous wave is wave 0 and has no effect.

  Psi1 <- Sigma - Delta_initial %*% Omega_full %*% t(Delta_initial)

  # Enforce positive semidefiniteness.

  Psi1 <- check_psd(Psi1)
  
  # Save Psi1.

  Psi_list[[1]] <- Psi1

  # M1 = Delta_initial %*% Omega, since the previous wave is wave 0 and has no effect.

  M1 <- Delta_initial %*% Omega_full

  # Save M1.

  M_list[[1]] <- M1

  # Simulate epsilon_1 with variance Psi1.

  e1 <- mvtnorm::rmvnorm(N, sigma = Psi1)

  # Simulate X1 and Y1 as C_full %*% t(Delta_initial) + e1.

  W1 <- C_full %*% t(Delta_initial) + e1

  # Save X1 and Y1.

  W_list[[1]] <- W1

  # ---- waves 2..T --------------------------------------------------------------------------------

  # Start with wave 2.

  for (t in 2:T_total) {

    # Extract Delta_t.

    Delta_t <- Delta_list[[t]]

    # Extract M from the previous wave.

    M_prev <- M_list[[t-1]]

    # Psi_t = Sigma - Phi %*% Sigma %*% t(Phi) -
    #         Delta_t %*% Omega %*% t(Delta_t) -
    #         Phi %*% M_prev %*% t(Delta_t) -
    #         Delta_t %*% t(M_prev) %*% t(Phi)

    Psi_t <-
      Sigma -
      Phi %*% Sigma %*% t(Phi) -
      Delta_t %*% Omega_full %*% t(Delta_t) -
      Phi %*% M_prev %*% t(Delta_t) -
      Delta_t %*% t(M_prev) %*% t(Phi)

    # Enforce positive semidefiniteness.

    Psi_t <- check_psd(Psi_t)

    # Save Psi.

    Psi_list[[t]] <- Psi_t

    # Update M:
    # M at wave t = Phi %*% M at wave t-1 + Delta_t %*% Omega.

    M_t <- Phi %*% M_prev + Delta_t %*% Omega_full

    # Save M.

    M_list[[t]] <- M_t

    # Simulate epsilon_t with variance Psi_t.

    e_t <- mvtnorm::rmvnorm(N, sigma = Psi_t)

    # Simulate X_t and Y_t as W_t = W_{t-1} %*% t(Phi) + C_full %*% t(Delta_t) + e_t.

    W_t <- W_list[[t-1]] %*% t(Phi) +
           C_full %*% t(Delta_t) +
           e_t

    # Save X_t and Y_t.

    W_list[[t]] <- W_t
  }

  # ---- build output data frame -------------------------------------------------------------------

  # Keep only the final T_obs waves after burn-in.

  keep_idx <- seq.int(from = burn_in + 1L, to = T_total)

  df <- matrix(NA, N, 2 * T_obs + ncol(C_full))

  colnames(df) <- c(
    paste0("x", 1:T_obs),
    paste0("y", 1:T_obs),
    feature_names
  )

  for (j in seq_along(keep_idx)) {
    t_keep <- keep_idx[j]
    df[, paste0("x", j)] <- W_list[[t_keep]][, 1]
    df[, paste0("y", j)] <- W_list[[t_keep]][, 2]
  }

  df[, (2 * T_obs + 1):(2 * T_obs + ncol(C_full))] <- C_full

  # ---- return object -----------------------------------------------------------------------------

  # The current engine returns only the final analysis data frame.
  # Any internal diagnostics used while constructing that data set stay internal.

  return(as.data.frame(df))
}