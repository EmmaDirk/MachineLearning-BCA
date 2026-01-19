# this function simulates panel data under a CLPM data-generating process, where  
# 1) we simulate the confounders first for each wave, based on the given Psi covariance matrix
# 2) we compute the variance these confounders induce at each wave given the B matrix for that wave
# 3) we compute how much variance must then be induced by the dynamic process to reach the target covariance 
# 4) however, we want some extra covariance between X and Y at each occasion, meaning that we need to first compute
#    the implied covariance at each wave, and then add this extra covariance to the target covariance at each wave
# 5) calculate the innovations covariance needed to reach this target covariance at each wave
# 6) simulate the panel data panel data for each wave using these innovations
#    then adding the lagged effects from A and the direct confounder effects from B
# 7) repeat for T waves with varying B matrices
#
# this script extends the simulate_panel_data function from the previous study by adding interaction effects
# between the confounders. Specifically, it generates interaction terms from the linear confounders,
# and adds them to the confounders. The interaction terms are orthogonalised with respect to their
# linear parents to ensure that they capture unique variance beyond the linear effects.
# ----------------------------------------------------------------------------------------------------------

simulate_panel_data_int <- function(
    N,                                                         # number of individuals
    T,                                                         # number of waves
    A,                                                         # 2x2 autoregressive/cross-lag matrix
    B_list,                                                    # list of B matrices: B_list[[t]] is 2 x p
    Psi,                                                       # k x k confounder covariance matrix (linear confounders)
    rho_extra,                                                 # extra covariance to add at observed level
    seed = NULL                                                # optional random seed for reproducibility
){

  # helper function to find stationary covariance c given A and S_U
  # given that the innovations are uncorrelated, what is the covariance between X and Y?
  # we need this to be able to add some extra covariance at the observed level on top of the existing covariance
  find_c <- function(A, S_U) {
    
    # given A and S_U, find the stationary covariance c between X and Y
    f <- function(c) {
      # given a candidate covariance c, compute the correlation of the innovations
      # we want this correlation to be 0, so we try to find a c that achieves this

      # target stationary covariance for (X_t, Y_t)
      # where the variance is 1 (as always) and covariance is c
      S_target_c <- matrix(c(1, c,
                             c, 1),
                           nrow = 2, byrow = TRUE)

      # dynamic component variance: S_dyn = S_target - S_U
      S_dyn_c <- S_target_c - S_U

      # ensure symmetry
      S_dyn_c <- (S_dyn_c + t(S_dyn_c)) / 2

      # innovations covariance implied by stationarity:
      # S_dyn = A S_dyn A' + Sigma_e_c
      Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A

      # ensure symmetry
      Sigma_e_c <- (Sigma_e_c + t(Sigma_e_c)) / 2

      # compute correlation of innovations (should match rho = 0)
      v1    <- Sigma_e_c[1, 1]
      v2    <- Sigma_e_c[2, 2]
      cov12 <- Sigma_e_c[1, 2]
      corr_e <- cov12 / sqrt(v1 * v2)

      corr_e - 0   # we want corr(e_x, e_y) = 0
    }

    # root-finding for covariance c between -0.99 and 0.99
    # wrapped in tryCatch so the simulation does not break if no root exists
    out <- tryCatch(
      uniroot(f, interval = c(-0.99, 0.99))$root,
      error = function(e) NA_real_
    )

    out
  }

  # set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # number of confounders (linear)
  k <- ncol(Psi)

  # simulate the confounders from their multivariate normal distribution
  # given a specified covariance matrix Psi, where diag=1 and for now the off-diag=0, giving I, the identity matrix. 
  U <- mvtnorm::rmvnorm(
    n     = N,                                                 # number of individuals
    mean  = rep(0, k),                                         # all confounders have mean 0
    sigma = Psi                                                # covariance matrix of confounders
  )

  # generate all interaction index sets of order >= 2
  int_idx <- if (k >= 2) {
    unlist(
      # for each order m from 2 to k
      # combine k confounders taken m at a time
      lapply(2:k, function(m) combn(1:k, m, simplify = FALSE)),
      # flatten the list
      recursive = FALSE
    )
  } else {
    list()
  }

  # build raw interaction terms
  U_int_raw <- if (length(int_idx) > 0) {
    sapply(int_idx, function(ix) {
      # compute product of confounders in this interaction
      apply(U[, ix, drop = FALSE], 1, prod)
    })
  } else {
    matrix(nrow = N, ncol = 0)
  }

  # orthogonalise interaction terms with respect to their linear parents
  # and normalise them to variance 1
  U_int_orth <- if (length(int_idx) > 0) {
    sapply(seq_along(int_idx), function(j) {

      # extract the parents for this interaction
      parents <- U[, int_idx[[j]], drop = FALSE]

      # raw interaction term
      z_raw <- U_int_raw[, j]

      # residualise the interaction on its parents
      z_orth <- resid(lm(z_raw ~ parents))

      # normalise to variance 1
      sdz <- sd(z_orth)
      if (is.na(sdz) || sdz == 0) z_orth else z_orth / sdz
    })
  } else {
    matrix(nrow = N, ncol = 0)
  }

  # combine linear and non-linear confounders
  U_full <- cbind(U, U_int_orth)

  
  # preparing containers for the variance structure
  S_dyn_list    <- vector("list", T)                           # variance coming from crosslaggs, autoreg and innovations
  Sigma_e_list  <- vector("list", T)                           # innovations covariance
  S_target_list <- vector("list", T)                           # target covariance at observed level
  c_base_vec    <- numeric(T)                                  # the baseline covariance if residual covariance is 0 (meaning all coming from system + confounders)
  c_total_vec   <- numeric(T)                                  # the total covariance at observed level (including extra rho)

  # computing the variance structure at each wave
  # for each wave t from 1 to T
  for (t in 1:T) {

    # get the B matrix for this wave
    B_t <- B_list[[t]]

    # the confounder induced variance covariance at this wave is B_t Psi B_t'
    # note: Psi is implicitly the identity for the full confounder set due to orthogonalisation
    S_U_t <- B_t %*% t(B_t)

    # the base covariance can be found (if the innovations are uncorrelated) by calling find_c
    c_base_t <- find_c(A, S_U_t)

    # store the base covariance in the container
    c_base_vec[t] <- c_base_t

    # now we need to add some extra covariance at the observed level
    c_total_t <- c_base_t + rho_extra

    # and store this in the container
    c_total_vec[t] <- c_total_t

    # S_target can then be specified using our computed total covariance and variance = 1
    S_target_t <- matrix(c(1, c_total_t,
                           c_total_t, 1),
                         nrow = 2, byrow = TRUE)
    
    # and save this in the container
    S_target_list[[t]] <- S_target_t

    # the variance coming from the dynamic process (cross-lag + auto + innovations) is then:
    # target - confounder induced variance
    S_dyn_t <- S_target_t - S_U_t

    # ensure symmetry
    S_dyn_t <- (S_dyn_t + t(S_dyn_t))/2   

    # store in container
    S_dyn_list[[t]] <- S_dyn_t

    # innovations covariance from stationarity:
    # S_dyn = A S_dyn A' + Sigma_e
    Sigma_e_t <- S_dyn_t - t(A) %*% S_dyn_t %*% A

    # ensure symmetry
    Sigma_e_t <- (Sigma_e_t + t(Sigma_e_t))/2

    # store in container
    Sigma_e_list[[t]] <- Sigma_e_t
  }

  # prepare data frame to hold simulated data
  # number of confounders
  p <- ncol(U_full)

  # data frame with N rows and 2*T + p columns
  df <- matrix(NA, nrow = N, ncol = 2*T + p)

  # set column names
  colnames(df) <- c(paste0("x", 1:T),
                    paste0("y", 1:T),
                    paste0("c", 1:p))

  # add confounders to dataframe
  df[, (2*T + 1):(2*T + p)] <- U_full

  # simulate the first wave, which is different because there are no lagged values yet
  Ddyn <- mvtnorm::rmvnorm(
    n     = N,                                                   # number of individuals
    mean  = c(0, 0),                                             # mean 0 for X and Y
    sigma = S_dyn_list[[1]]                                      # variance covariance matrix at wave 1
  ) 

  # add the direct confounder effects
  obs1 <- Ddyn + U_full %*% t(B_list[[1]])

  # store in dataframe
  df[, "x1"] <- obs1[, 1]
  df[, "y1"] <- obs1[, 2]

  # simulate waves 2 to T
  for (t in 2:T) {

    # pull variance covariance matrix for this wave
    Sigma_e_t <- Sigma_e_list[[t]]

    # dynamic process
    Ddyn <- Ddyn %*% t(A) + mvtnorm::rmvnorm(N, sigma = Sigma_e_t)

    # add direct confounder effects
    obs <- Ddyn + U_full %*% t(B_list[[t]])

    # store
    df[, paste0("x", t)] <- obs[, 1]
    df[, paste0("y", t)] <- obs[, 2]
  }

  return(df)
}