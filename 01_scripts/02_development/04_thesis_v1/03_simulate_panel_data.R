# This function simulates panel data under a cross-lagged data-generating process
# with time-invariant confounders whose effects may vary over time.
#
# The goal is to generate two observed variables (X and Y) measured over T waves
# for N individuals, where X and Y are confounded by a time-invariant set of
# confounders C (including interaction features).
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
# Omega11    = covariance matrix of base confounders (c1..ck)
# Delta_list = list of confounder effect matrices Delta_t
#
# Delta_t maps confounder features into (X_t, Y_t):
#
#   W_t = Phi W_{t-1} + Delta_t C + epsilon_t
#
# where
#
#   W_t = (X_t, Y_t)'
#   C   = vector of confounder features
#   epsilon_t ~ (0, Psi_t)
#
# Because the same confounders affect every wave,
# the lagged term Phi W_{t-1} and the direct confounder term Delta_t C
# are correlated. Therefore we track:
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
# ------------------------------------------------------------------------------------

simulate_panel_data <- function(
    N,                                                            # number of individuals
    T,                                                            # number of waves
    Phi,                                                          # lag matrix (diagonal elements are autoregressive effects)
    Delta_list,                                                   # list of Delta matrices
    Omega11,                                                      # covariance matrix of base confounders
    Sigma,                                                        # desired covariance matrix of (X_t, Y_t)
    seed = NULL,                                                  # seed for reproducibility
    eig_tol = 1e-10                                               # tolerance for positive semidefiniteness
){

  # ------------------------- input checks -------------------------

  # check that N and T are positive integers
  if (!is.numeric(N) || N < 1 || N != as.integer(N))
    stop("N must be a positive integer.")
  if (!is.numeric(T) || T < 1 || T != as.integer(T))
    stop("T must be a positive integer.")

  # check that Phi and Sigma are 2x2 matrices
  if (!is.matrix(Phi) || !all(dim(Phi) == c(2,2)))
    stop("Phi must be a 2x2 matrix.")
  if (!is.matrix(Sigma) || !all(dim(Sigma) == c(2,2)))
    stop("Sigma must be a 2x2 matrix.")

  # check that sigma is symmetric
  if (!isTRUE(all.equal(Sigma, t(Sigma))))
    stop("Sigma must be symmetric.")
  
  # check that sigma has 1 on the diagonal
  if (!isTRUE(all.equal(diag(Sigma), c(1,1))))
    stop("Sigma must have 1 on the diagonal.")

  # check that Omega11 is a square matrix
  if (!is.matrix(Omega11))
    stop("Omega11 must be a matrix.")

  # check that Omega11 is symmetric
  if (!isTRUE(all.equal(Omega11, t(Omega11))))
    stop("Omega11 must be symmetric.")

  # check that Omega11 has 1 on the diagonal
  if (!isTRUE(all.equal(diag(Omega11), rep(1, nrow(Omega11)))))
    stop("Base confounders must be standardized (diag(Omega11)=1).")

  # check that Delta_list is a list
  if (!is.list(Delta_list) || length(Delta_list) != T)
    stop("Delta_list must be a list of length T.")

  # check that Delta matrices have column names
  feature_names <- colnames(Delta_list[[1]])
  if (is.null(feature_names))
    stop("Delta matrices must have column names.")

  # k = number of base confounders
  k <- nrow(Omega11)

  # with linear names of the form c1:cK
  lin_names <- paste0("c",1:k)

  # ------------------------- helper: PSD check -------------------------
  # checks that covariance matrices are positive semidefinite

  check_psd <- function(S) {

    # take mean of transpose and itself
    S <- (S + t(S))/2

    # eigen decomposition
    evd <- eigen(S, symmetric = TRUE)

    # extract eigenvalues
    vals <- evd$values

    # if any eigenvalue < eig_tol, stop
    if (any(vals < -eig_tol))
      stop("A covariance matrix is not positive semidefinite.")

    # otherwise, return the original symmetrized matrix
    S
  }

  # ------------------------- helper: parse feature name -------------------------
  # converts "c1:2:3" -> c(1,2,3)

  parse_feature <- function(name) {

    ix <- as.integer(strsplit(sub("^c","",name),":")[[1]])
    ix
  }

  # ------------------------- helper: analytic covariance of confounder features -------------------------
  #
  # We do not need (but sure can) to estimate Cov(C_full) empirically from the simulated sample.
  # Instead, because all features are either:
  # - standardized main effects
  # - standardized two-way interactions
  # - standardized three-way interactions
  #
  # and because the base confounders are multivariate normal with covariance Omega11,
  # we can derive the full covariance matrix analytically.
  #
  # This gives the population covariance of the feature vector C,
  # which is more stable than cov(C_full) and does not depend on N.

  # for centered Gaussian variables, the sixth moment equals the sum over all 15 pairings
  # of the 6 indices. we build those 15 pairings once and reuse them.
  pairings6 <- local({

    rec <- function(v) {

      # if there are no indices left, return one empty pairing structure
      if (length(v) == 0) return(list(list()))

      # if there are exactly two indices left, they form the final pair
      if (length(v) == 2) return(list(list(c(v[1], v[2]))))

      # otherwise recursively pair the first element with each later element
      first <- v[1]
      out <- list()

      for (m in 2:length(v)) {

        # remove the chosen pair and recurse on the remainder
        rest <- v[-c(1, m)]
        sub  <- rec(rest)

        # append the chosen pair to each recursive solution
        for (s in sub) {
          out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
        }
      }

      out
    }

    rec(1:6)
  })

  # compute E[X1 X2 X3 X4 X5 X6] for centered Gaussian variables
  # using Isserlis' theorem: sum over all pairings of products of covariances
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

  # variance of the raw two-way product C_i * C_l
  # for standardized Gaussian variables this is 1 + rho^2
  var_raw_2way <- function(i, l, Omega11) {

    1 + Omega11[i,l]^2
  }

  # variance of the raw three-way product C_i * C_l * C_m
  # this is the denominator used to standardize the three-way interaction
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

  # covariance between two main effects
  cov_main_main <- function(a, b, Omega11) {

    Omega11[a,b]
  }

  # covariance between a main effect and a standardized two-way interaction
  # for centered Gaussian variables this is 0
  cov_main_2way <- function(a, i, l, Omega11) {

    0
  }

  # covariance between a main effect and a standardized three-way interaction
  cov_main_3way <- function(a, i, l, m, Omega11) {

    num <- Omega11[a,i] * Omega11[l,m] +
           Omega11[a,l] * Omega11[i,m] +
           Omega11[a,m] * Omega11[i,l]

    den <- sqrt(var_raw_3way(i, l, m, Omega11))

    num / den
  }

  # covariance between two standardized two-way interactions
  cov_2way_2way <- function(i, l, p, q, Omega11) {

    num <- Omega11[i,p] * Omega11[l,q] +
           Omega11[i,q] * Omega11[l,p]

    den <- sqrt(
      var_raw_2way(i, l, Omega11) *
      var_raw_2way(p, q, Omega11)
    )

    num / den
  }

  # covariance between a standardized two-way interaction and a standardized three-way interaction
  # for centered Gaussian variables this is 0
  cov_2way_3way <- function(i, l, p, q, r, Omega11) {

    0
  }

  # covariance between two standardized three-way interactions
  cov_3way_3way <- function(i, l, m, p, q, r, Omega11) {

    num <- sixth_moment_gaussian(c(i, l, m, p, q, r), Omega11)

    den <- sqrt(
      var_raw_3way(i, l, m, Omega11) *
      var_raw_3way(p, q, r, Omega11)
    )

    num / den
  }

  # determine what kind of feature we are dealing with
  # check that the integer is between 1 and 6
  feature_type <- function(ix) {

    order <- length(ix)

    if (order == 1) return("main")
    if (order == 2) return("2way")
    if (order == 3) return("3way")

    stop("Only up to three-way interactions supported.")
  }

  # compute the covariance between any two supported features
  # this works directly with parsed feature indices
  cov_feature_pair <- function(ix1, ix2, Omega11) {

    type1 <- feature_type(ix1)
    type2 <- feature_type(ix2)

    # main with main
    if (type1 == "main" && type2 == "main") {
      return(cov_main_main(ix1[1], ix2[1], Omega11))
    }

    # main with 2way
    if (type1 == "main" && type2 == "2way") {
      return(cov_main_2way(ix1[1], ix2[1], ix2[2], Omega11))
    }

    # 2way with main
    if (type1 == "2way" && type2 == "main") {
      return(cov_main_2way(ix2[1], ix1[1], ix1[2], Omega11))
    }

    # main with 3way
    if (type1 == "main" && type2 == "3way") {
      return(cov_main_3way(ix1[1], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # 3way with main
    if (type1 == "3way" && type2 == "main") {
      return(cov_main_3way(ix2[1], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # 2way with 2way
    if (type1 == "2way" && type2 == "2way") {
      return(cov_2way_2way(ix1[1], ix1[2], ix2[1], ix2[2], Omega11))
    }

    # 2way with 3way
    if (type1 == "2way" && type2 == "3way") {
      return(cov_2way_3way(ix1[1], ix1[2], ix2[1], ix2[2], ix2[3], Omega11))
    }

    # 3way with 2way
    if (type1 == "3way" && type2 == "2way") {
      return(cov_2way_3way(ix2[1], ix2[2], ix1[1], ix1[2], ix1[3], Omega11))
    }

    # 3way with 3way
    if (type1 == "3way" && type2 == "3way") {
      return(cov_3way_3way(ix1[1], ix1[2], ix1[3], ix2[1], ix2[2], ix2[3], Omega11))
    }

    stop("Unsupported feature combination.")
  }

  # build the full covariance matrix of the confounder feature vector
  # in exactly the same order as feature_names
  build_full_Omega_from_features <- function(feature_names, Omega11) {

    # p = total number of requested features
    p <- length(feature_names)

    # parse every feature name once
    parsed_features <- lapply(feature_names, parse_feature)

    # initialize Omega
    Omega <- matrix(0, p, p)
    rownames(Omega) <- feature_names
    colnames(Omega) <- feature_names

    # fill the upper triangle and mirror it to the lower triangle
    for (a in seq_len(p)) {
      for (b in a:p) {

        val <- cov_feature_pair(parsed_features[[a]], parsed_features[[b]], Omega11)

        Omega[a,b] <- val
        Omega[b,a] <- val
      }
    }

    Omega
  }

  # ------------------------- helper: construct confounder feature matrix -------------------------
  #
  # Base confounders are simulated directly.
  # Interaction terms are constructed analytically and standardized using
  # the same formulas assumed in sample_delta_1().

  build_confounder_features <- function(C_base, feature_names, Omega11) {

    # N = number  base confounders
    N <- nrow(C_base)

    # initialize the matrix of confounders
    out <- matrix(NA, N, length(feature_names))
    colnames(out) <- feature_names

    # for every entry in the feature matrix
    for (j in seq_along(feature_names)) {

      # get the name of the feature
      name <- feature_names[j]

      # parse to get the interaction order
      ix <- parse_feature(name)

      # the order is now the number of elements
      order <- length(ix)

      # main effect
      if (order == 1) {

        # if order is 1, simply return the base confounder
        out[,j] <- C_base[,ix]

      # two-way interaction
      } else if (order == 2) {

        # make the index
        i <- ix[1]; l <- ix[2]

        # extract the interaction between the base confounders
        # E[C_iC_l]
        rho <- Omega11[i,l]

        # compute the raw interaction
        raw <- C_base[,i] * C_base[,l]

        # now standardize using the formula
        out[,j] <- (raw - rho) / sqrt(1 + rho^2)

      # three-way interaction
      } else if (order == 3) {

        # make the index
        i <- ix[1]; l <- ix[2]; m <- ix[3]

        # get the covariances: E[C_iC_l], E[C_iC_m], E[C_lC_m]
        rho_il <- Omega11[i,l]
        rho_im <- Omega11[i,m]
        rho_lm <- Omega11[l,m]

        # make the raw interaction
        raw <- C_base[,i] * C_base[,l] * C_base[,m]

        # compute the variance of that raw interaction
        denom <- sqrt(
          1 +
          2*rho_il^2 +
          2*rho_im^2 +
          2*rho_lm^2 +
          8*rho_il*rho_im*rho_lm
        )

        # standardize the output
        out[,j] <- raw / denom

      } else {

        stop("Only up to three-way interactions supported.")
      }
    }

    out
  }

  # ------------------------- seed -------------------------

  if (!is.null(seed))
    set.seed(seed)

  # ------------------------- simulate base confounders -------------------------

  C_base <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(0,k),
    sigma = Omega11
  )

  colnames(C_base) <- lin_names

  # ------------------------- build full confounder feature matrix -------------------------

  C_full <- build_confounder_features(
    C_base,
    feature_names,
    Omega11
  )

  # analytic covariance of the confounder feature matrix
  Omega_full <- build_full_Omega_from_features(
    feature_names,
    Omega11
  )

  # ------------------------- containers -------------------------

  Psi_list <- vector("list",T)
  M_list <- vector("list",T)
  W_list <- vector("list",T)

  # ------------------------- wave 1 -------------------------

  # extract Delta1
  Delta1 <- Delta_list[[1]]

  # Psi1 = Sigma - Delta1 %*% Omega %*% t(Delta1)
  # since the previous wave is wave 0 and has no effect
  Psi1 <- Sigma - Delta1 %*% Omega_full %*% t(Delta1)

  # enforce positive semidefiniteness
  Psi1 <- check_psd(Psi1)
  
  # save psi1
  Psi_list[[1]] <- Psi1

  # M1 = Delta1 %*% Omega
  # since the previous wave is wave 0 and has no effect
  M1 <- Delta1 %*% Omega_full

  # save M1
  M_list[[1]] <- M1

  # now we can simulate epsilon 1, with variance Psi1
  e1 <- mvtnorm::rmvnorm(N, sigma = Psi1)

  # we can simulate X1 and Y1 as: C_full %*% t(Delta1) + e1
  W1 <- C_full %*% t(Delta1) + e1

  # save X1 and Y1
  W_list[[1]] <- W1

  # ------------------------- waves 2..T -------------------------

  # start with wave 2
  for (t in 2:T) {

    # extract Delta_t
    Delta_t <- Delta_list[[t]]

    # extract M from previous wave
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

    # enforce positive semidefiniteness
    Psi_t <- check_psd(Psi_t)

    # save Psi
    Psi_list[[t]] <- Psi_t

    # update M
    # M at wave t = Phi %*% M at wave t-1 + Delta_t %*% Omega
    M_t <- Phi %*% M_prev + Delta_t %*% Omega_full

    # save M
    M_list[[t]] <- M_t

    # simulate epsilon_t with variance Psi_t
    e_t <- mvtnorm::rmvnorm(N, sigma = Psi_t)

    # simulate X_t and Y_t as W_t = W_{t-1} %*% t(Phi) + C_full %*% t(Delta_t) + e_t
    W_t <- W_list[[t-1]] %*% t(Phi) +
           C_full %*% t(Delta_t) +
           e_t

    # save X_t and Y_t
    W_list[[t]] <- W_t
  }

  # ------------------------- build output data frame -------------------------

  df <- matrix(NA, N, 2*T + ncol(C_full))

  colnames(df) <- c(
    paste0("x",1:T),
    paste0("y",1:T),
    feature_names
  )

  for (t in seq_len(T)) {

    df[,paste0("x",t)] <- W_list[[t]][,1]
    df[,paste0("y",t)] <- W_list[[t]][,2]
  }

  df[,(2*T+1):(2*T+ncol(C_full))] <- C_full

  # ------------------------- return object -------------------------

  return(list(

    data = as.data.frame(df),

    confounders = as.data.frame(C_full),

    base_confounders = as.data.frame(C_base),

    Psi_list = Psi_list,

    M_list = M_list,

    Sigma = Sigma,

    Phi = Phi,

    Omega11 = Omega11,

    Omega_full = Omega_full,

    Delta_list = Delta_list
  ))
}