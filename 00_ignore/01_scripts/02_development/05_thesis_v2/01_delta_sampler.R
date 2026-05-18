# -----------------------------------------------------------------------------------------------------------------
# This script contains a function to sample the coefficients with which time-invariant confounders (e.g. age),
# affect the variables X and Y at time t = 1.
# Since the autoregessive effects are denoted using beta, and the cross-lagged effects are denoted using gamma,
# we denote the confounder effects using delta. As such, here we sample the delta coefficients at time t = 1.
# do note that t=1 is the first time point in the DGM, but with burn-in, this is not the first time point
# that is analysed. 
#
# The function samples the coefficients such that:
# - the total variance explained by the confounders at time t = 1 is R2_total
# - the fraction of variance explained by interaction terms is R2_interaction
#
# The argument k refers to the number of base confounders (main effects): c1..ck.
# Interaction terms are only included when requested via toggles:
# - include_2way: include every possible two-way interaction (e.g. c1:c2, c1:c3, ...)
# - include_3way: include every possible three-way interaction (e.g. c1:c2:c4, ...)
#
# This is accomplished by sampling random coefficients from a normal distribution, and scaling them
# such that the variance explained by the linear and interaction parts equals the desired R^2 values.
# The function checks whether the absolute values of the sampled coefficients are within the desired range.
# If not, it resamples until valid coefficients are found, or the maximum number of tries is reached.
#
# This requires the following assumptions / restictions:
# - All variables and X and Y have an expected value of 0.
# - All variables have been standardized to have a variance of 1.
# - All variables are normally distributed.
# -----------------------------------------------------------------------------------------------------------------

sample_delta_1 <- function(
  k,                             # number of base confounders
  Omega11,                       # k x k covariance matrix of the base confounders
  R2_total,                      # target Var(eta); if Var(X)=1, this equals the target R^2
  R2_interaction = 0,            # fraction of R2_total allocated to the nonlinear bucket
  include_2way = FALSE,          # include all two-way interaction terms
  include_3way = FALSE,          # include all three-way interaction terms
  force_positive = FALSE,        # ensure all sampled coefficients are positive
  min_abs = 0,                   # minimum allowed absolute coefficient size
  max_abs = 1,                   # maximum allowed absolute coefficient size
  max_tries = 100000,            # maximum number of resampling attempts
  eig_tol = 1e-10                # tolerance for positive semidefiniteness check
) {

  # ------------------------- input checks -------------------------
  #
  # This part of the script simply checks whether the input arguments are valid.

  # if Omega11 is not a k x k matrix, throw an error
  if (!is.matrix(Omega11) || nrow(Omega11) != k || ncol(Omega11) != k)
    stop("Omega11 must be a k x k matrix.")

  # if Omega11 is not symmetric, throw an error
  if (!isTRUE(all.equal(Omega11, t(Omega11), tolerance = 1e-12)))
    stop("Omega11 must be symmetric.")

  # if Omega11 is not positive semidefinite, throw an error
  ev <- eigen(Omega11, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < -eig_tol))
    stop("Omega11 must be positive semidefinite.")

  # if Omega11 does not have ones on the diagonal, throw an error
  # the base confounders are assumed to be standardized
  if (!isTRUE(all.equal(diag(Omega11), rep(1, k), tolerance = 1e-12)))
    stop("Omega11 must have 1 on the diagonal because the base confounders are standardized.")

  # if R2_total is not between 0 and 1, throw an error
  if (R2_total < 0 || R2_total > 1)
    stop("R2_total must be between 0 and 1.")

  # if R2_interaction is not between 0 and 1, throw an error
  if (R2_interaction < 0 || R2_interaction > 1)
    stop("R2_interaction must be between 0 and 1.")

  # if two-way interactions are requested, we need at least 2 base confounders
  if (include_2way && k < 2)
    stop("include_2way requires k >= 2.")

  # if three-way interactions are requested, we need at least 3 base confounders
  if (include_3way && k < 3)
    stop("include_3way requires k >= 3.")

  # if interaction terms are requested, interaction variance must be > 0
  if ((include_2way || include_3way) && R2_interaction == 0)
    stop("Interaction variance must be > 0 when interaction terms are included.")

  # if interaction variance is > 0, at least one interaction toggle must be TRUE
  if (R2_interaction > 0 && !(include_2way || include_3way))
    stop("Interaction variance is positive but no interaction terms are included.")

  # if min_abs is negative, throw an error
  if (min_abs < 0)
    stop("min_abs must be >= 0.")

  # if max_abs is less than min_abs, throw an error
  if (max_abs < min_abs)
    stop("max_abs must be >= min_abs.")

  # ------------------------- dimensions ---------------------------
  #
  # This part of the script computes the number of possible interaction terms.

  # if we include interaction terms, we need to compute the number of possible interaction terms.
  # if we include 2way interaction, its number is k choose 2
  m2 <- if (include_2way) choose(k, 2) else 0

  # if we include 3way interaction, its number is k choose 3
  m3 <- if (include_3way) choose(k, 3) else 0

  # then the total number of interaction terms is m2 + m3
  m_int <- m2 + m3

  # targets
  # The target variance of the linear component is (1 - R2_interaction) * R2_total
  Vlin_star <- (1 - R2_interaction) * R2_total

  # The target variance of the interaction component is R2_interaction * R2_total
  VNL_star  <- R2_interaction * R2_total

  # helpers
  # makes the index names: 1:2, 1:3, 2:3, ...
  all_pairs <- function(k) t(combn(seq_len(k), 2))

  # makes the index names: 1:2:3, ...
  all_triples <- function(k) t(combn(seq_len(k), 3))

  # ------------------------- covariance objects -------------------------
  #
  # In this part we compute every part of the covariance matrix Omega. 

  # These functions work directly with the final standardized variables:
  # - Ci       = standardized base confounder i
  # - Z_ij     = standardized two-way interaction between Ci and Cj
  # - T_ijl    = standardized three-way interaction between Ci, Cj, and Cl
  #
  # it should be noted that:
  # - Omega11:  Cov(Ci, Cj) = Omega11_ij
  # - Omega12:  Cov(Ci, Z_ij) = 0
  # - Omega22:  Cov(Z_ij, Z_pq)
  # - Omega13:  Cov(Ci, T_ijl)
  # - Omega33:  Cov(T_ijl, T_pqr)

  # Helper: 
  # for centered Gaussian variables, the sixth moment equals the sum over all 15 pairings.
  # we precompute those 15 pairings once and reuse them.
  # this is basically just a helper to build the indices
  pairings6 <- local({

    # for every possible pair
    rec <- function(v) {

      # if there are no pairs left, return an empty list
      if (length(v) == 0) return(list(list()))
      
      # if there is only one pair left, return a list with that pair
      if (length(v) == 2) return(list(list(c(v[1], v[2]))))

      # otherwise, recursively build the list
      # pick the first element
      first <- v[1]

      # initialize the list
      out <- list()

      # pair the first element with the rest
      for (m in 2:length(v)) {

        # remove the first and the current element
        rest <- v[-c(1, m)]

        # do this recursively
        sub  <- rec(rest)

        # for each sub-list
        for (s in sub) {

          # add the first and the current element
          out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
        }
      }
      # return the list
      out
    }
    # call the recursive function
    rec(1:6)
  })

  # Helper:
  # compute E[X1 X2 X3 X4 X5 X6] for centered Gaussian variables
  # takes a vector of indices and a covariance matrix of the base confounders
  sixth_moment_gaussian <- function(idx, Omega11) {

    # loop over 15 pairings
    # take the sum of the covariance between the two indices
    sum(vapply(
      pairings6,

      function(pr) {

        # find the covariance between the two indices
        prod(vapply(
          pr,
          function(pa) Omega11[idx[pa[1]], idx[pa[2]]],

          # where each iteration returns a number
          numeric(1)
        ))
      },
      numeric(1)
    ))
  }

  # analytic covariance between two standardized two-way interactions
  # this is written directly as Cov(Z_ij, Z_pq)
  # where i and j are the indices of the first interaction
  # and p and q are the indices of the second interaction
  cov_Z2 <- function(i, j, p, q, Omega11) {

    # the numerator is the sum of the product of the covariances
    num <- Omega11[i, p] * Omega11[j, q] +
           Omega11[i, q] * Omega11[j, p]

    # the denominator is the standardizing factor
    den <- sqrt((1 + Omega11[i, j]^2) *
                (1 + Omega11[p, q]^2))

    # the covariance is the numerator divided by the denominator
    num / den
  }

  # analytic covariance between two standardized three-way interactions
  # this is written directly as Cov(T_ijl, T_pqr)
  # it takes six indices, three for each interaction
  cov_T3 <- function(i, j, l, p, q, r, Omega11) {

    # the numerator is the sum of the product of the covariances
    # which we can calculate analytically using the sixth moment function
    # as: E[C1 C2 C3 C4 C5 C6]
    num <- sixth_moment_gaussian(c(i, j, l, p, q, r), Omega11)

    # the denominator is the standardizing factor
    # which is just the variance of the two interactions
    # which can be computed analytically:
    den1 <- 1 +
      2 * Omega11[i, j]^2 +
      2 * Omega11[i, l]^2 +
      2 * Omega11[j, l]^2 +
      8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]

    den2 <- 1 +
      2 * Omega11[p, q]^2 +
      2 * Omega11[p, r]^2 +
      2 * Omega11[q, r]^2 +
      8 * Omega11[p, q] * Omega11[p, r] * Omega11[q, r]

    num / sqrt(den1 * den2)
  }

  # analytic covariance between a standardized main effect and a standardized three-way interaction
  # this is written directly as Cov(C_a, T_ijl)
  # takes the index of the main effect and the indices of the three-way interaction
  # and the covariance matrix of the base confounders
  cov_13 <- function(a, i, j, l, Omega11) {

    # the numerator is the sum of the product of the covariances
    num <- Omega11[a, i] * Omega11[j, l] +
           Omega11[a, j] * Omega11[i, l] +
           Omega11[a, l] * Omega11[i, j]

    # the denominator is the standardizing factor
    den <- sqrt(
      1 +
      2 * Omega11[i, j]^2 +
      2 * Omega11[i, l]^2 +
      2 * Omega11[j, l]^2 +
      8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]
    )

    num / den
  }

  # we start by building the covariance matrix of the two-way interactions
  build_Omega22 <- function(Omega11) {

    # get all possible pairs
    P <- all_pairs(nrow(Omega11))

    # they are m2 pairs
    m2 <- nrow(P)

    # initialize the m2 x m2 matrix Omega
    Omega22 <- matrix(0, m2, m2)

    # if there are no pairs, return Omega empty
    if (m2 == 0) return(Omega22)

    # since all two-way interaction variables are standardized, the diagonal is 1
    diag(Omega22) <- 1

    # for every pair
    for (a in seq_len(m2)) {

      # get the indices to select the a'th pair
      i <- P[a, 1]; j <- P[a, 2]

      # for every other pair, compute the covariance analytically
      for (b in seq_len(m2)) {

        # the diagonal is already known
        if (a == b) next

        # get the indices to select the b'th pair
        p <- P[b, 1]; q <- P[b, 2]

        # fill the off-diagonal entry directly as Cov(Z_ij, Z_pq)
        Omega22[a, b] <- cov_Z2(i, j, p, q, Omega11)
      }
    }

    Omega22
  }

  # build the covariance matrix of the three-way interactions
  build_Omega33 <- function(Omega11) {

    # get all possible triples
    T <- all_triples(nrow(Omega11))

    # they are m3 triples
    m3 <- nrow(T)

    # initialize the m3 x m3 matrix Omega
    Omega33 <- matrix(0, m3, m3)

    # if there are no triples, return Omega empty
    if (m3 == 0) return(Omega33)

    # since all three-way interaction variables are standardized, the diagonal is 1
    diag(Omega33) <- 1

    # for every triple
    for (a in seq_len(m3)) {

      # get the indices to select the a'th triple
      i <- T[a, 1]; j <- T[a, 2]; l <- T[a, 3]

      # for every other triple, compute the covariance analytically
      for (b in seq_len(m3)) {

        # the diagonal is already known
        if (a == b) next

        # get the indices to select the b'th triple
        p <- T[b, 1]; q <- T[b, 2]; r <- T[b, 3]

        # fill the off-diagonal entry directly as Cov(T_ijl, T_pqr)
        Omega33[a, b] <- cov_T3(i, j, l, p, q, r, Omega11)
      }
    }

    Omega33
  }

  # build the covariance matrix between the main effects and the three-way interactions
  build_Omega13 <- function(Omega11) {

    # get all possible triples
    T <- all_triples(nrow(Omega11))

    # they are m3 triples
    m3 <- nrow(T)

    # initialize the k x m3 matrix Omega
    Omega13 <- matrix(0, k, m3)

    # if there are no triples, return Omega empty
    if (m3 == 0) return(Omega13)

    # for every triple
    for (t in seq_len(m3)) {

      # get the indices to select the t'th triple
      i <- T[t, 1]; j <- T[t, 2]; l <- T[t, 3]

      # for every main effect, compute the covariance analytically
      for (a in seq_len(k)) {

        # fill the entry directly as Cov(C_a, T_ijl)
        Omega13[a, t] <- cov_13(a, i, j, l, Omega11)
      }
    }

    Omega13
  }

  # ---------------------------- solve scales -------------------------------
  # this function finds two scaling parameters sL (for main effects) and s (for interactions)
  # to transform the sampled coefficients to the desired R^2 values. 
  # A = t(bL) %*% Omega11 %*% bL
  # B = t(b2) %*% Omega22 %*% b2
  # C = t(b3) %*% Omega33 %*% b3
  # D = t(bL) %*% Omega13 %*% b3
  # which are computed before calling this function

  solve_scales <- function(A, B, C, D, Vlin_star, VNL_star) {

    # if target linear variance is zero, all linear coefficients must be zero
    if (Vlin_star == 0) {

      # do note that then the interaction variance must be > 0. 
      if ((B + C) <= 0)
        stop("B + C must be > 0 when interaction variance is positive.")
      
      # then sL = 0, and s = sqrt(VNL_star / (B + C))
      return(list(sL = 0, s = sqrt(VNL_star / (B + C)), r = Inf))
    }

    # if target interaction variance is zero, all interaction coefficients must be zero
    if (VNL_star == 0) {

      # then there must be linear variance
      if (A <= 0)
        stop("A must be > 0.")
      
      # then sL = sqrt(Vlin_star / A), and s = 0
      return(list(sL = sqrt(Vlin_star / A), s = 0, r = 0))
    }

    # similarly, if there are no interaction terms
    if ((B + C) == 0 && D == 0) {

      # then A must be > 0
      if (A <= 0) stop("A must be > 0.")
      
      # and sL = sqrt(Vlin_star / A), and s = 0
      return(list(sL = sqrt(Vlin_star / A), s = 0))
    }

    # if there is no covariance between the linear and 3-way interaction
    if (D == 0) {

      # then A must be > 0
      if (A <= 0) stop("A must be > 0.")
      
      # and B + C must be > 0
      if ((B + C) <= 0) stop("B + C must be > 0.")
      
      # then sL = sqrt(Vlin_star / A), and s = sqrt(VNL_star / (B + C))
      return(list(
        sL = sqrt(Vlin_star / A),
        s  = sqrt(VNL_star / (B + C))
      ))
    }

    # otherwise, solve the quadratic equation (abc-formula)
    # note down each part of the formula
    a <- Vlin_star * (B + C)
    b <- D * (Vlin_star - VNL_star)
    c <- -VNL_star * A

    # compute the discriminant
    disc <- b^2 - 4 * a * c

    # check that the discriminant is positive
    if (disc < 0)
      stop("No real solution for the ratio r = s / sL.")

    # compute the two roots
    r1 <- (-b + sqrt(disc)) / (2 * a)
    r2 <- (-b - sqrt(disc)) / (2 * a)

    # check that the roots are admissible
    admissible <- function(r) is.finite(r) && (r > 0) && (A + r * D > 0)

    # choose the admissible root
    r <- if (admissible(r1)) r1 else if (admissible(r2)) r2 else NA_real_

    # again check that the root is admissible
    if (!is.finite(r))
      stop("No admissible positive solution for r = s / sL.")

    # and use the admissible root to compute sL and s
    sL <- sqrt(Vlin_star / (A + r * D))
    s  <- r * sL

    list(sL = sL, s = s, r = r)
  }

  #  ----------------------- covariance blocks ----------------------------
  # these depend only on Omega11, so we compute them once.
  # initialize the Omega matrices
  Omega22 <- matrix(0, m2, m2)
  Omega33 <- matrix(0, m3, m3)
  Omega13 <- matrix(0, k,  m3)

  # compute the Omega matrices if toggled
  if (include_2way) Omega22 <- build_Omega22(Omega11)
  if (include_3way) {
    Omega33 <- build_Omega33(Omega11)
    Omega13 <- build_Omega13(Omega11)
  }

  # names 
  # if toggled, we need to name the Omega matrices
  # function to generate the names for the 2-way interactions
  int2_names <- if (include_2way) {
    combn(1:k, 2, FUN = function(ix) paste0("c", paste(ix, collapse = ":")))
  } else character(0)

  # function to generate the names for the 3-way interactions
  int3_names <- if (include_3way) {
    combn(1:k, 3, FUN = function(ix) paste0("c", paste(ix, collapse = ":")))
  } else character(0)

  # names for the linear terms
  lin_names <- paste0("c", 1:k)

  # all feature names
  feature_names <- c(lin_names, int2_names, int3_names)

  # if toggled, name the Omega matrices
  if (m2 > 0) {
    rownames(Omega22) <- int2_names
    colnames(Omega22) <- int2_names
  }

  if (m3 > 0) {
    rownames(Omega33) <- int3_names
    colnames(Omega33) <- int3_names
    rownames(Omega13) <- lin_names
    colnames(Omega13) <- int3_names
  }

  # full Omega 
  # now we can put it all together
  # this function takes Omega11, Omega22, Omega33, and Omega13
  build_full_Omega <- function(Omega11, Omega22, Omega33, Omega13) {

    # compute the total number of features
    p_total <- k + nrow(Omega22) + nrow(Omega33)

    # initialize Omega
    Omega <- matrix(0, p_total, p_total)

    # split Omega into blocks
    idx_L  <- seq_len(k)
    idx_2  <- if (m2 > 0) (k + 1):(k + m2) else integer(0)
    idx_3  <- if (m3 > 0) (k + m2 + 1):(k + m2 + m3) else integer(0)

    # Omega11 is the covariance matrix of the standardized main effects
    Omega[idx_L, idx_L] <- Omega11

    # the covariance matrix of the standardized two-way interaction terms
    if (m2 > 0) Omega[idx_2, idx_2] <- Omega22

    # the covariance matrix of the standardized three-way interaction terms
    if (m3 > 0) Omega[idx_3, idx_3] <- Omega33

    # the covariance matrix between the standardized main effects and standardized three-way interactions
    if (m3 > 0) {
      Omega[idx_L, idx_3] <- Omega13
      Omega[idx_3, idx_L] <- t(Omega13)
    }

    # it should be noted that:
    # - the covariance of main effects with two-way interactions is zero
    # - the covariance of two-way interactions with three-way interactions is zero

    rownames(Omega) <- feature_names
    colnames(Omega) <- feature_names
    Omega
  }

  # ------------------------- sample one row -------------------------
  sample_one_row <- function() {

    # sample a random direction for the linear coefficients
    bL <- rnorm(k)

    # sample a random direction for the interaction coefficients
    b_int <- if (m_int > 0) rnorm(m_int) else numeric(0)

    # if toggled, force all sampled coefficients to be positive
    if (force_positive) {
      bL <- abs(bL)
      if (m_int > 0) b_int <- abs(b_int)
    }

    # split the interaction coefficients into the 2way and 3way parts
    b2 <- if (m2 > 0) b_int[1:m2] else numeric(0)
    b3 <- if (m3 > 0) b_int[(m2 + 1):(m2 + m3)] else numeric(0)

    # A = variance contribution of the linear part
    A <- as.numeric(t(bL) %*% Omega11 %*% bL)

    # B = variance contribution of the 2way interaction part
    B <- if (m2 > 0) {
      as.numeric(t(b2) %*% Omega22 %*% b2)
    } else 0

    # C = variance contribution of the 3way interaction part
    C <- if (m3 > 0) {
      as.numeric(t(b3) %*% Omega33 %*% b3)
    } else 0

    # D = covariance contribution between the linear and 3way parts
    D <- if (m3 > 0) {
      as.numeric(t(bL) %*% Omega13 %*% b3)
    } else 0

    # if there are no interaction terms, scale only the linear part
    if (m_int == 0) {
      if (A <= 0) stop("A must be > 0; check Omega11.")
      sL <- sqrt(Vlin_star / A)
      s  <- 0
    } else {
      scales <- solve_scales(
        A = A, B = B, C = C, D = D,
        Vlin_star = Vlin_star,
        VNL_star  = VNL_star
      )
      sL <- scales$sL
      s  <- scales$s
    }

    # scale the sampled directions into the final coefficient vectors
    dL <- sL * bL
    d2 <- if (m2 > 0) s * b2 else numeric(0)
    d3 <- if (m3 > 0) s * b3 else numeric(0)

    # return the full sampled coefficient vector
    c(dL, d2, d3)
  }

  # ------------------------- main sampling loop -------------------------
  for (i in seq_len(max_tries)) {

    dx <- sample_one_row()
    dy <- sample_one_row()

    if (all(abs(c(dx, dy)) >= min_abs & abs(c(dx, dy)) <= max_abs)) {

      D1 <- rbind(dx, dy)
      rownames(D1) <- c("X", "Y")
      colnames(D1) <- feature_names

      Omega_full <- build_full_Omega(Omega11, Omega22, Omega33, Omega13)

      return(list(
        Delta = D1,
        Omega = Omega_full,
        Omega_blocks = list(
          Omega11 = Omega11,
          Omega22 = Omega22,
          Omega33 = Omega33,
          Omega13 = Omega13
        )
      ))
    }
  }

  stop("Failed to sample Delta at t=1 within max_tries. Try relaxing bounds or increasing max_tries.")
}