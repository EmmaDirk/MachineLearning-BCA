# this script contains functions to generate delta trajectories from our baseline sampled Delta-Matrix D1
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
generate_D_constant <- function(
    D1,                                                           # baseline D-matrix
    T                                                             # number of time points
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(D1))
    stop("D1 must be a matrix.")

  if (T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # create an emtpy list of length T
  D_list <- vector("list", T)

  # set names for each time point
  names(D_list) <- paste0("t", 1:T)

  # fill the list with copies of D1
  for (t in 1:T) {
    D_list[[t]] <- D1
  }
  
  # return the list of D matrices
  return(D_list)
}


# function (2) to generate D trajectory with a hard step
generate_D_stepwise <- function(
    D1,                                                            # baseline D-matrix
    T,                                                             # number of time points
    step_at = floor(T/2) + 1,                                      # when the step starts (default: second half)
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40                                                  # higher (or lower) R2 (after the step)
){

  # ------------------------- input checks -------------------------
  if (!is.matrix(D1))
    stop("D1 must be a matrix.")

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

  # build list of D matrices of length T
  D_list <- vector("list", T)

  # set names for each time point
  names(D_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      D_list[[t]] <- D1

    } else {

      # if we are AT the step or AFTER the step: jump to the higher delta matrix
      # this creates the hard step like:
      # {D1, D1, D1, D1*scale_factor, D1*scale_factor}
      D_list[[t]] <- D1 * scale_factor
    }
  }

  # return the list of D matrices
  return(D_list)
}

# function (3) to generate D trajectory with a hard step,
# where the post-step coefficients are a convex mixture of the baseline pattern
# and a random alternative pattern, followed by exact rescaling so that
# the variance decomposition into linear and nonlinear parts is preserved
generate_D_stepwise_mixture <- function(
    D1,                                                            # baseline D-matrix
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
  if (!is.matrix(D1))
    stop("D1 must be a matrix.")

  # check that Omega is a matrix
  if (!is.matrix(Omega))
    stop("Omega must be a matrix.")

  # check that D1 and Omega are conformable
  if (ncol(D1) != nrow(Omega) || ncol(D1) != ncol(Omega))
    stop("D1 and Omega are not conformable: ncol(D1) must equal nrow(Omega) = ncol(Omega).")

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
  if (is.null(colnames(D1)))
    stop("D1 must have column names so linear and nonlinear terms can be identified.")

  # optional reproducibility
  if (!is.null(seed)) set.seed(seed)

  # ------------------------- identify coefficient blocks -------------------------
  # separate the linear and nonlinear coefficients
  # helper to count how many separators (colons) are in a string
  count_colons <- function(x) {
    nchar(gsub("[^:]", "", x))
  }

  # the number of colons in each coefficient name
  n_colons <- count_colons(colnames(D1))

  # the indices of the linear and nonlinear coefficients
  idx_L  <- which(n_colons == 0)                                   # c1, c2, c3, ...
  idx_NL <- which(n_colons >= 1)                                   # c1:c2, c1:c2:c3, ...

  # check that D1 has linear coefficients
  if (length(idx_L) == 0)
    stop("No linear coefficients found in D1.")

  # do note: nonlinear coefficients are optional here
  has_NL <- length(idx_NL) > 0

  # linear and nonlinear covariance matrices
  Omega_L  <- Omega[idx_L, idx_L, drop = FALSE]

  # where the non-linear coefficients are present, otherwise NULL
  Omega_NL <- if (has_NL) Omega[idx_NL, idx_NL, drop = FALSE] else NULL

  # ------------------------- helper to compute quadratic-form variance -------------------------
  # compute the variance of a vector deltas (d) under a covariance matrix (Om)
  qvar <- function(d, Om) {

    # as: d' Om d
    as.numeric(t(d) %*% Om %*% d)
  }

  # ------------------------- helper to build one post-step row -------------------------
  make_post_row <- function(d_old_row) {

    # split old row into linear and nonlinear parts
    dL_old  <- d_old_row[idx_L]
    dNL_old <- if (has_NL) d_old_row[idx_NL] else numeric(0)

    # compute baseline variance decomposition
    V_L_old  <- qvar(dL_old, Omega_L)
    V_NL_old <- if (has_NL) qvar(dNL_old, Omega_NL) else 0
    V_old    <- V_L_old + V_NL_old

    # check that old_R2 matches the baseline row variance implied by D1
    if (abs(V_old - old_R2) > 1e-8)
      stop("The variance implied by D1 does not match old_R2 for at least one row.")

    # compute target linear and nonlinear variances after the step
    if (has_NL) {
      share_L  <- V_L_old / V_old
      share_NL <- V_NL_old / V_old

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
    scale_L  <- sqrt(V_L_new / V_L_mix)
    scale_NL <- if (has_NL) {
      if (V_NL_new <= tol) 0 else sqrt(V_NL_new / V_NL_mix)
    } else {
      0
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
  D_post <- D1

  for (r in seq_len(nrow(D1))) {
    D_post[r, ] <- make_post_row(D1[r, ])
  }

  rownames(D_post) <- rownames(D1)
  colnames(D_post) <- colnames(D1)

  # ------------------------- build list of D matrices of length T -------------------------
  D_list <- vector("list", T)

  # set names for each time point
  names(D_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      D_list[[t]] <- D1

    } else {

      # if we are AT the step or AFTER the step: jump to the mixed-and-rescaled delta matrix
      D_list[[t]] <- D_post
    }
  }

  # return the list of D matrices
  return(D_list)
}
