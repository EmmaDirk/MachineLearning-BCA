# this script contains a function to sample beta coefficients for confounders
# such that the total R^2 of confounders at time t = 1 is equal to a specified value
# and such that a fraction eta_1 of that R^2 is allocated to non-linear (interaction) terms
# ---------------------------------------------------------------------

sample_B1 <- function(
  k,                                                             # number of confounders
  R2_1,                                                          # total confounder R^2 at t = 1
  eta_1     = 0.0,                                               # fraction of R2_1 allocated to non-linear terms
  min_abs   = 0.01,                                              # minimum absolute value for each beta
  max_abs   = 0.60,                                              # maximum absolute value for each beta
  max_tries = 100000                                             # maximum sampling attempts
) {

  # check eta_1 is a valid fraction
  if (eta_1 < 0 || eta_1 > 1)
    stop("eta_1 must be between 0 and 1")

  # number of linear terms (always k main effects)
  k_lin <- k

  # number of non-linear terms (all interaction terms of order 2..k)
  k_non <- if (eta_1 > 0) {
    sum(sapply(2:k, function(m) choose(k, m)))
  } else 0

  # split total R^2 into linear and non-linear parts
  R2_lin <- (1 - eta_1) * R2_1
  R2_non <- eta_1 * R2_1

  # start the loop: for i in the max number of tries
  for (i in seq_len(max_tries)) {

    # ----------------------------
    # sample LINEAR effects (X)
    # ----------------------------

    # if we have linear terms and nonzero linear R^2, sample + normalize to hit R2_lin
    bx_lin <- if (k_lin > 0 && R2_lin > 0) {

      # sample k_lin random numbers from a normal distribution
      v <- rnorm(k_lin)

      # normalize so the sum of squares of this vector is 1,
      # then scale so the sum of squares is R2_lin
      sqrt(R2_lin) * v / sqrt(sum(v^2))

    } else numeric(0)  # otherwise, no linear effects

    # ----------------------------
    # sample LINEAR effects (Y)
    # ----------------------------

    by_lin <- if (k_lin > 0 && R2_lin > 0) {

      # sample k_lin random numbers from a normal distribution
      v <- rnorm(k_lin)

      # normalize + scale to hit R2_lin
      sqrt(R2_lin) * v / sqrt(sum(v^2))

    } else numeric(0)

    # initialize non-linear vectors as empty
    bx_non <- by_non <- numeric(0)

    # ----------------------------
    # sample NON-LINEAR effects (X, Y)
    # ----------------------------

    if (k_non > 0 && R2_non > 0) {

      # sample k_non random numbers from a normal distribution for X
      v <- rnorm(k_non)

      # normalize + scale to hit R2_non
      bx_non <- sqrt(R2_non) * v / sqrt(sum(v^2))

      # sample k_non random numbers from a normal distribution for Y
      v <- rnorm(k_non)

      # normalize + scale to hit R2_non
      by_non <- sqrt(R2_non) * v / sqrt(sum(v^2))
    }

    # combine linear + non-linear parts into full beta vectors
    bx <- c(bx_lin, bx_non)
    by <- c(by_lin, by_non)

    # check if all absolute values are within the specified bounds
    if (all(abs(c(bx, by)) >= min_abs & abs(c(bx, by)) <= max_abs)) {

      # if so, return the B matrix (2 rows: X and Y)
      B1 <- rbind(bx, by)

      # set row names
      rownames(B1) <- c("X", "Y")

      # set column names:
      # - first k are main effects: c1..ck
      # - then (if present) interaction names like c1:2, c1:3, ..., c1:2:3, ...
      if (k_non == 0) {

        # only linear terms
        colnames(B1) <- paste0("c", 1:k)

      } else {

        # build interaction names in the same order as combn(1:k, m)
        int_names <- unlist(
          lapply(2:k, function(m)
            combn(1:k, m, FUN = function(ix)
              paste0("c", paste(ix, collapse=":")))
          )
        )

        # concatenate main effect names + interaction names
        colnames(B1) <- c(paste0("c", 1:k), int_names)
      }

      return(B1)
    }
  }

  # if the loop doesn't return a valid B matrix, throw an error
  stop("Failed to sample B1 within max_tries")
}
