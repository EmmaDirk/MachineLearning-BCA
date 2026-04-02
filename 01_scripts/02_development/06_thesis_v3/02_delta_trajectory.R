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
#
# IMPORTANT:
# T is now interpreted as the number of observed waves the user wants to analyse.
# burn_in is added automatically inside the trajectory generators.
# So if step_at = 3 and burn_in = 20, the step is placed at internal wave 23.
# ---------------------------------------------------------------------

# ------------------------------------ helper ------------------------------------
# validate the shared trajectory inputs
validate_Delta_trajectory_inputs <- function(T, burn_in) {

  # check that T is a positive integer
  if (T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # check that burn_in is a non-negative integer
  if (burn_in < 0 || burn_in != as.integer(burn_in))
    stop("burn_in must be a non-negative integer.")
}


# helper to prepend burn-in copies of Delta1 to the observed trajectory
prepend_burn_in <- function(Delta_obs_list, Delta1, burn_in) {

  # if there is no burn-in, just return the observed trajectory
  if (burn_in == 0L) {
    return(Delta_obs_list)
  }

  # copy the original sampled coefficient matrix through the burn-in
  burn_list <- vector("list", burn_in)

  # fill the burn-in list with copies of Delta1
  for (b in seq_len(burn_in)) {
    burn_list[[b]] <- Delta1
  }

  # combine burn-in and observed trajectory
  c(burn_list, Delta_obs_list)
}


# ------------------------------------ 1) constant effects ------------------------------------
generate_Delta_constant <- function(
    Delta1,                                                       # baseline Delta-matrix
    T,                                                            # number of observed time points
    burn_in = 0L                                                  # number of burn-in waves
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(Delta1))
    stop("Delta1 must be a matrix.")

  validate_Delta_trajectory_inputs(T = T, burn_in = burn_in)

  # coerce to integers
  T <- as.integer(T)
  burn_in <- as.integer(burn_in)

  # create an emtpy list of length T
  Delta_obs_list <- vector("list", T)

  # fill the observed-wave list with copies of Delta1
  for (t in 1:T) {
    Delta_obs_list[[t]] <- Delta1
  }

  # prepend burn-in copies automatically
  Delta_list <- prepend_burn_in(
    Delta_obs_list = Delta_obs_list,
    Delta1 = Delta1,
    burn_in = burn_in
  )

  # set names for each internal time point
  names(Delta_list) <- paste0("t", seq_along(Delta_list))
  
  # return the list of Delta matrices
  return(Delta_list)
}


# ------------------------------------ 2) stepwise effects --------------------------------------
generate_Delta_stepwise <- function(
    Delta1,                                                        # baseline Delta-matrix
    T,                                                             # number of observed time points
    burn_in = 0L,                                                  # number of burn-in waves
    step_at = floor(T/2) + 1,                                      # when the step starts in observed time
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40                                                  # higher (or lower) R2 (after the step)
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(Delta1))
    stop("Delta1 must be a matrix.")

  validate_Delta_trajectory_inputs(T = T, burn_in = burn_in)

  # coerce to integers
  T <- as.integer(T)
  burn_in <- as.integer(burn_in)

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
  Delta_obs_list <- vector("list", T)

  # fill the observed-wave list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      Delta_obs_list[[t]] <- Delta1

    } else {

      # if we are AT the step or AFTER the step: jump to the higher delta matrix
      # this creates the hard step like:
      # {D1, D1, D1, D1*scale_factor, D1*scale_factor}
      Delta_obs_list[[t]] <- Delta1 * scale_factor
    }
  }

  # prepend burn-in copies automatically
  Delta_list <- prepend_burn_in(
    Delta_obs_list = Delta_obs_list,
    Delta1 = Delta1,
    burn_in = burn_in
  )

  # set names for each internal time point
  names(Delta_list) <- paste0("t", seq_along(Delta_list))

  # return the list of Delta matrices
  return(Delta_list)
}

# ------------------------------------ 3) stepwise mixture effects --------------------------------------
# where the post-step coefficients are a convex mixture of the baseline pattern
# and a random alternative pattern, followed by exact rescaling so that
# the variance decomposition into linear and nonlinear parts is preserved
generate_Delta_stepwise_mixture <- function(
    Delta1,                                                        # baseline Delta-matrix
    T,                                                             # number of observed time points
    Omega,                                                         # full covariance matrix of all confounder features
    burn_in = 0L,                                                  # number of burn-in waves
    step_at = floor(T/2) + 1,                                      # when the step starts in observed time
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

  # check that T and burn_in are valid
  validate_Delta_trajectory_inputs(T = T, burn_in = burn_in)

  # coerce to integers
  T <- as.integer(T)
  burn_in <- as.integer(burn_in)

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

  # helper
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

    # otherwise, solve the quadratic equation in r = sNL / sL
    # derived from:
    #   V_L_star / V_NL_star = (A + rD) / (r^2 B + rD)
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

    # check which root is admissible
    admissible <- function(r) is.finite(r) && (r > 0) && ((A + r * D) > tol)

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

    # split the old row into linear and nonlinear parts
    dL_old  <- d_old_row[idx_L]
    dNL_old <- if (has_NL) d_old_row[idx_NL] else numeric(0)

    # compute old variance decomposition
    V_L_old  <- qvar(dL_old, Omega_L)
    V_NL_old <- if (has_NL) qvar(dNL_old, Omega_NL) else 0
    V_cross  <- if (has_NL) as.numeric(t(dL_old) %*% Omega_LNL %*% dNL_old) else 0

    # split the cross term equally between the two blocks
    V_L_eff_old  <- V_L_old + V_cross
    V_NL_eff_old <- V_NL_old + V_cross

    # total explained variance should equal old_R2
    V_old_total <- V_L_eff_old + V_NL_eff_old - V_cross

    # check this up to numerical tolerance
    if (abs(V_old_total - old_R2) > 1e-6)
      stop("The supplied Delta1 row does not match old_R2 under the supplied Omega.")

    # preserve the old relative variance split, but rescale total variance to new_R2
    if (old_R2 > tol) {
      share_L  <- V_L_eff_old  / old_R2
      share_NL <- V_NL_eff_old / old_R2
    } else {
      share_L  <- 1
      share_NL <- 0
    }

    V_L_star  <- share_L  * new_R2
    V_NL_star <- share_NL * new_R2

    # sample random alternative directions
    dL_new_raw  <- rnorm(length(dL_old))
    dNL_new_raw <- if (has_NL) rnorm(length(dNL_old)) else numeric(0)

    # build convex mixtures
    dL_mix  <- (1 - lambda_L)  * dL_old  + lambda_L  * dL_new_raw
    dNL_mix <- if (has_NL) (1 - lambda_NL) * dNL_old + lambda_NL * dNL_new_raw else numeric(0)

    # compute the variance pieces for the mixed directions
    A <- qvar(dL_mix, Omega_L)
    B <- if (has_NL) qvar(dNL_mix, Omega_NL) else 0
    D <- if (has_NL) as.numeric(t(dL_mix) %*% Omega_LNL %*% dNL_mix) else 0

    # solve for the rescaling constants
    if (has_NL) {
      scales <- solve_block_scales(
        A = A,
        B = B,
        D = D,
        V_L_star = V_L_star,
        V_NL_star = V_NL_star
      )
      dL_post  <- scales$sL  * dL_mix
      dNL_post <- scales$sNL * dNL_mix
    } else {
      # linear-only case
      if (A <= tol)
        stop("A must be > 0 in the linear-only case.")
      dL_post  <- sqrt(new_R2 / A) * dL_mix
      dNL_post <- numeric(0)
    }

    # reconstruct the full coefficient row
    out <- numeric(length(d_old_row))
    out[idx_L] <- dL_post
    if (has_NL) out[idx_NL] <- dNL_post
    out
  }

  # ------------------------- create the post-step matrix -------------------------
  Delta_post <- Delta1

  for (r in seq_len(nrow(Delta1))) {
    Delta_post[r, ] <- make_post_row(Delta1[r, ])
  }

  # build list of Delta matrices of length T
  Delta_obs_list <- vector("list", T)

  # fill the observed-wave list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      Delta_obs_list[[t]] <- Delta1

    } else {

      # at and after the step, use the mixed-and-rescaled post-step matrix
      Delta_obs_list[[t]] <- Delta_post
    }
  }

  # prepend burn-in copies automatically
  Delta_list <- prepend_burn_in(
    Delta_obs_list = Delta_obs_list,
    Delta1 = Delta1,
    burn_in = burn_in
  )

  # set names for each internal time point
  names(Delta_list) <- paste0("t", seq_along(Delta_list))

  # return the list of Delta matrices
  return(Delta_list)
}