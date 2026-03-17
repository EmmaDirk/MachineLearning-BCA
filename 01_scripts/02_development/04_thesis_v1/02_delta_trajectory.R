# this script contains functions to generate delta trajectories from our baseline sampled Delta-Matrix Delta1
# these delta trajectories represent how the effects of baseline confounders change over time
# we have two scenarios implemented here:
# - the effects remain constant over time (do not change)
# - the effects change in a stepwise manner at a specified time point
# - the effects change in a stepwise manner at a specified time point, but their rank order may change
#
# function 1: the effects of baseline confounders remain constant over time
# function 2: the effects of baseline confounders change in a stepwise manner over time such that 
#             the variance explained (R2) by the confounders increases from old_R2 to new_R2 at a 
#             specified time point
# function 3: expands on function 2 by allowing the rank order of the confounders to change. We do this
#             by sampling a random vector of coefficients. Then we generate new confounder coefficients as:
#             delta_bar = (1-lambda)*delta_bar_old + lambda*delta_bar_new
#             when labda = 1, the old confounder efffects are fully replaced by the new confounder effects.
#             finally, we need to rescale the coefficients such that they guarantee new R2 = new_R2. 
# ---------------------------------------------------------------------

# function (1) to generate Delta trajectory where the coefficients are constant
generate_Delta_constant <- function(
    Delta1,                                                       # baseline Delta-matrix
    T                                                             # number of time points
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(Delta1))
    stop("Delta1 must be a matrix.")

  if (T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # create an emtpy list of length T
  Delta_list <- vector("list", T)

  # set names for each time point
  names(Delta_list) <- paste0("t", 1:T)

  # fill the list with copies of Delta1
  for (t in 1:T) {
    Delta_list[[t]] <- Delta1
  }
  
  # return the list of Delta matrices
  return(Delta_list)
}


# function (2) to generate Delta trajectory with a hard step
generate_Delta_stepwise <- function(
    Delta1,                                                        # baseline Delta-matrix
    T,                                                             # number of time points
    step_at = floor(T/2) + 1,                                      # when the step starts (default: second half)
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40                                                  # higher (or lower) R2 (after the step)
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(Delta1))
    stop("Delta1 must be a matrix.")

  if (T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  if (step_at < 1 || step_at > T || step_at != as.integer(step_at))
    stop("step_at must be an integer between 1 and T.")

  if (!is.numeric(old_R2) || length(old_R2) != 1 || old_R2 <= 0 || old_R2 > 1)
    stop("old_R2 must be a single number in (0, 1].")

  if (!is.numeric(new_R2) || length(new_R2) != 1 || new_R2 < 0 || new_R2 > 1)
    stop("new_R2 must be a single number in [0, 1].")

  # scaling factor such that R2 changes from old_R2 to new_R2
  # scaling factor = sqrt(new_R2 / old_R2)
  scale_factor <- sqrt(new_R2 / old_R2)                            # for old_R2=0.15 and new_R2=0.40 -> ~ 1.633

  # build list of Delta matrices of length T
  Delta_list <- vector("list", T)

  # set names for each time point
  names(Delta_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      Delta_list[[t]] <- Delta1

    } else {

      # if we are AT the step or AFTER the step: jump to the higher delta matrix
      # this creates the hard step like:
      # {D1, D1, D1, D1*scale_factor, D1*scale_factor}
      Delta_list[[t]] <- Delta1 * scale_factor
    }
  }

  # return the list of Delta matrices
  return(Delta_list)
}
# function (3) to generate Delta trajectory with a hard step,
# where the post-step coefficients are a convex mixture of the baseline pattern
# and a random alternative pattern, followed by exact rescaling so that
# the variance decomposition into linear and nonlinear parts is preserved
generate_Delta_stepwise_mixture <- function(
    Delta1,                                                        # baseline Delta-matrix
    T,                                                             # number of time points
    Omega,                                                         # full covariance matrix of all confounder features
    step_at = floor(T/2) + 1,                                      # when the step starts (default: second half)
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40,                                                 # higher (or lower) R2 (after the step)
    lambda_L = 0.50,                                               # mixture weight for linear coefficients
    lambda_NL = 0.50,                                              # mixture weight for nonlinear coefficients
    seed = NULL,                                                   # optional seed for reproducibility
    tol = 1e-10                                                    # numerical tolerance
){

  # ------------------------- input checks -------------------------
  # check that D1 is a matrix
  if (!is.matrix(Delta1))
    stop("Delta1 must be a matrix.")

  # check that Omega is a matrix
  if (!is.matrix(Omega))
    stop("Omega must be a matrix.")

  # check that D1 and Omega are conformable
  if (ncol(Delta1) != nrow(Omega) || ncol(Delta1) != ncol(Omega))
    stop("Delta1 and Omega are not conformable: ncol(Delta1) must equal nrow(Omega) = ncol(Omega).")

  # check that Omega is symmetric
  if (!isTRUE(all.equal(Omega, t(Omega), tolerance = tol)))
    stop("Omega must be symmetric.")

  # check that T is a positive integer
  if (T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # check that step_at is an integer between 1 and T
  if (step_at < 1 || step_at > T || step_at != as.integer(step_at))
    stop("step_at must be an integer between 1 and T.")

  # check that old_R2 and new_R2 are numbers in (0, 1]
  if (!is.numeric(old_R2) || length(old_R2) != 1 || old_R2 <= 0 || old_R2 > 1)
    stop("old_R2 must be a single number in (0, 1].")

  # check that old_R2 and new_R2 are numbers in (0, 1]
  if (!is.numeric(new_R2) || length(new_R2) != 1 || new_R2 < 0 || new_R2 > 1)
    stop("new_R2 must be a single number in [0, 1].")

  # check that lambda_L and lambda_NL are numbers in [0, 1]
  if (!is.numeric(lambda_L) || length(lambda_L) != 1 || lambda_L < 0 || lambda_L > 1)
    stop("lambda_L must be a single number in [0, 1].")

  # check that lambda_L and lambda_NL are numbers in [0, 1]
  if (!is.numeric(lambda_NL) || length(lambda_NL) != 1 || lambda_NL < 0 || lambda_NL > 1)
    stop("lambda_NL must be a single number in [0, 1].")

  # check that D1 has column names
  if (is.null(colnames(Delta1)))
    stop("Delta1 must have column names so linear and nonlinear terms can be identified.")

  # optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  # ------------------------- identify coefficient blocks -------------------------
  # separate the linear and nonlinear coefficients
  # helper to count how many separators (colons) are in a string
  count_colons <- function(x) {
    nchar(gsub("[^:]", "", x))
  }

  # the number of colons in each coefficient name
  n_colons <- count_colons(colnames(Delta1))

  # the indices of the linear and nonlinear coefficients
  idx_L  <- which(n_colons == 0)                                   # c1, c2, c3, ...
  idx_NL <- which(n_colons >= 1)                                   # c1:c2, c1:c2:c3, ...

  # check that D1 has linear coefficients
  if (length(idx_L) == 0)
    stop("No linear coefficients found in Delta1.")

  # do note: nonlinear coefficients are optional here
  has_NL <- length(idx_NL) > 0

  # linear and nonlinear covariance matrices
  Omega_L  <- Omega[idx_L, idx_L, drop = FALSE]

  # where the non-linear coefficients are present, otherwise NULL
  Omega_NL <- if (has_NL) Omega[idx_NL, idx_NL, drop = FALSE] else NULL

  # covariance matrix between the linear and nonlinear blocks
  Omega_LNL <- if (has_NL) Omega[idx_L, idx_NL, drop = FALSE] else NULL

  # ------------------------- helper to compute quadratic-form variance -------------------------
  # compute the variance of a vector deltas (d) under a covariance matrix (Om)
  qvar <- function(d, Om) {

    # as: d' Om d
    as.numeric(t(d) %*% Om %*% d)
  }

  # ------------------------- helper to solve block scales -------------------------
  # this function finds two scaling parameters sL (for linear effects) and sNL (for nonlinear effects)
  # such that:
  # - the linear effective variance hits V_L_star
  # - the nonlinear effective variance hits V_NL_star
  #
  # where the linear and nonlinear blocks each receive half of the linear-nonlinear cross covariance:
  # - linear effective variance     = sL^2 * A + sL * sNL * D
  # - nonlinear effective variance  = sNL^2 * B + sL * sNL * D
  #
  # A = dL'  Omega_L   dL
  # B = dNL' Omega_NL  dNL
  # D = dL'  Omega_LNL dNL
  solve_block_scales <- function(A, B, D, V_L_star, V_NL_star) {

    # if target linear variance is zero, all linear coefficients must be zero
    if (V_L_star <= tol) {

      # then the nonlinear block must absorb all remaining variance
      if (B <= tol)
        stop("B must be > 0 when nonlinear target variance is positive.")

      # then sL = 0, and sNL = sqrt(V_NL_star / B)
      return(list(sL = 0, sNL = sqrt(V_NL_star / B), r = Inf))
    }

    # if target nonlinear variance is zero, all nonlinear coefficients must be zero
    if (V_NL_star <= tol) {

      # then the linear block must absorb all remaining variance
      if (A <= tol)
        stop("A must be > 0 when linear target variance is positive.")

      # then sL = sqrt(V_L_star / A), and sNL = 0
      return(list(sL = sqrt(V_L_star / A), sNL = 0, r = 0))
    }

    # if there is no covariance between the linear and nonlinear blocks
    if (abs(D) <= tol) {

      # then A must be > 0
      if (A <= tol) stop("A must be > 0.")

      # and B must be > 0
      if (B <= tol) stop("B must be > 0.")

      # then the two blocks can be scaled separately
      return(list(
        sL  = sqrt(V_L_star  / A),
        sNL = sqrt(V_NL_star / B)
      ))
    }

    # otherwise, solve the quadratic equation (abc-formula)
    # note down each part of the formula
    a <- V_L_star * B
    b <- D * (V_L_star - V_NL_star)
    c <- -V_NL_star * A

    # compute the discriminant
    disc <- b^2 - 4 * a * c

    # check that the discriminant is positive
    if (disc < 0)
      stop("No real solution for the ratio r = sNL / sL.")

    # compute the two roots
    r1 <- (-b + sqrt(disc)) / (2 * a)
    r2 <- (-b - sqrt(disc)) / (2 * a)

    # check that the roots are admissible
    admissible <- function(r) {
      is.finite(r) && (r > 0) && ((A + r * D) > tol) && ((r^2 * B + r * D) > tol)
    }

    # choose the admissible root
    r <- if (admissible(r1)) r1 else if (admissible(r2)) r2 else NA_real_

    # again check that the root is admissible
    if (!is.finite(r))
      stop("No admissible positive solution for r = sNL / sL.")

    # and use the admissible root to compute sL and sNL
    sL  <- sqrt(V_L_star / (A + r * D))
    sNL <- r * sL

    list(sL = sL, sNL = sNL, r = r)
  }

  # ------------------------- helper to build one post-step row -------------------------
  make_post_row <- function(d_old_row) {

    # split old row into linear and nonlinear parts
    dL_old  <- d_old_row[idx_L]
    dNL_old <- if (has_NL) d_old_row[idx_NL] else numeric(0)

    # compute baseline variance decomposition
    V_L_old  <- qvar(dL_old, Omega_L)
    V_NL_old <- if (has_NL) qvar(dNL_old, Omega_NL) else 0
    V_cross_old_half <- if (has_NL) as.numeric(t(dL_old) %*% Omega_LNL %*% dNL_old) else 0
    V_L_eff_old  <- V_L_old  + V_cross_old_half
    V_NL_eff_old <- V_NL_old + V_cross_old_half
    V_old        <- V_L_eff_old + V_NL_eff_old

    # check that old_R2 matches the baseline row variance implied by D1
    if (abs(V_old - old_R2) > 1e-8)
      stop("The variance implied by Delta1 does not match old_R2 for at least one row.")

    # compute target linear and nonlinear variances after the step
    if (has_NL) {
      share_L  <- V_L_eff_old  / V_old
      share_NL <- V_NL_eff_old / V_old

      V_L_new  <- share_L  * new_R2
      V_NL_new <- share_NL * new_R2
    } else {
      V_L_new  <- new_R2
      V_NL_new <- 0
    }

    # sample random alternative directions
    uL  <- rnorm(length(dL_old))
    uNL <- if (has_NL) rnorm(length(dNL_old)) else numeric(0)

    # create convex mixtures
    dL_mix  <- (1 - lambda_L)  * dL_old  + lambda_L  * uL
    dNL_mix <- if (has_NL) (1 - lambda_NL) * dNL_old + lambda_NL * uNL else numeric(0)

    # check that mixed directions are not degenerate
    V_L_mix  <- qvar(dL_mix, Omega_L)
    V_NL_mix <- if (has_NL) qvar(dNL_mix, Omega_NL) else 0

    # if the mixed linear direction has near-zero variance, try again
    if (V_L_mix <= tol)
      stop("The mixed linear direction has near-zero variance. Try a different seed or lambda_L.")

    if (has_NL && V_NL_mix <= tol && V_NL_new > tol)
      stop("The mixed nonlinear direction has near-zero variance. Try a different seed or lambda_NL.")

    # rescale mixtures so that the target variance decomposition is exactly preserved
    if (has_NL) {

      # define the quadratic-form pieces for the mixed directions
      A <- qvar(dL_mix, Omega_L)
      B <- qvar(dNL_mix, Omega_NL)
      D <- as.numeric(t(dL_mix) %*% Omega_LNL %*% dNL_mix)

      # solve for the block-specific scales
      scales <- solve_block_scales(
        A = A,
        B = B,
        D = D,
        V_L_star = V_L_new,
        V_NL_star = V_NL_new
      )

      scale_L  <- scales$sL
      scale_NL <- scales$sNL

    } else {

      # if there are no nonlinear coefficients, only the linear block remains
      scale_L  <- sqrt(V_L_new / V_L_mix)
      scale_NL <- 0
    }

    dL_new  <- scale_L * dL_mix
    dNL_new <- if (has_NL) scale_NL * dNL_mix else numeric(0)

    # combine back into one full row
    d_new <- numeric(length(d_old_row))
    d_new[idx_L] <- dL_new

    if (has_NL) {
      d_new[idx_NL] <- dNL_new
    }

    d_new
  }

  # ------------------------- create the post-step D matrix -------------------------
  Delta_post <- Delta1

  for (r in seq_len(nrow(Delta1))) {
    Delta_post[r, ] <- make_post_row(Delta1[r, ])
  }

  rownames(Delta_post) <- rownames(Delta1)
  colnames(Delta_post) <- colnames(Delta1)

  # ------------------------- build list of Delta matrices of length T -------------------------
  Delta_list <- vector("list", T)

  # set names for each time point
  names(Delta_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      Delta_list[[t]] <- Delta1

    } else {

      # if we are AT the step or AFTER the step: jump to the mixed-and-rescaled delta matrix
      Delta_list[[t]] <- Delta_post
    }
  }

  # return the list of Delta matrices
  return(Delta_list)
}