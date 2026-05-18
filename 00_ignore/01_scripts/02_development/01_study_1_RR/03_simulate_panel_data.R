# This function simulates panel data under a CLPM (cross-lagged panel model) data-generating process:
# its goal is to produce data for two observed variables (X and Y) measured over T waves for N individuals,
# where X and Y are confounded by a time-invariant set of confounders U.
# It aims to achieve an observed variance of 1 for all variables, such that the effects can be interpreted as standardized effects.
#
# We take the following steps: 
# 1) simulate time-invariant confounders U for N individuals from a multivariate normal distribution with covariance matrix Psi
# 2) for each wave t, use the wave-specific delta matrix D_t to compute the confounder-induced covariance: S_U_t = D_t Psi D_t'
# 3) determine the stationary covariance between X and Y implied by the system (given A and S_U_t), assuming uncorrelated innovations
# 4) add extra observed-level covariance (rho_extra) to obtain the final target covariance of X and Y at wave t. This will show up as residual covariance in the CLPM. 
# 5) compute the remaining covariance that must come from the dynamic component: S_dyn_t = S_target_t - S_U_t
# 6) derive the innovation covariance Sigma_e_t that makes the process stationary at wave t
# 7) simulate (X_t, Y_t) over T waves by iterating the dynamic system via A and adding direct confounder effects via D_t
# ------------------------------------------------------------------------------------------------------------

simulate_panel_data <- function(
    N,                                                         # number of individuals
    T,                                                         # number of waves
    A,                                                         # 2x2 autoregressive/cross-lag matrix
    D_list,                                                    # list of D matrices: D_list[[t]] is 2 x k
    Psi,                                                       # k x k confounder covariance matrix
    rho_extra,                                                 # extra covariance to add at observed level. 
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
      Sigma_e_c <- S_dyn_c - A %*% S_dyn_c %*% t(A)

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

    # return the found covariance c
    out
  }

  # set seed, and use the supplied seed if there is one
  if (!is.null(seed)) set.seed(seed)

  # number of confounders
  k <- ncol(Psi)

  # simulate the confounders from their multivariate normal distribution
  # given a specified covariance matrix Psi, where diag=1 and for now the off-diag=0, giving I, the identity matrix. 
  U <- mvtnorm::rmvnorm(
    n     = N,                                                 # number of individuals
    mean  = rep(0, k),                                         # all confounders have mean 0
    sigma = Psi                                                # covariance matrix of confounders
  )

  # preparing containers for the variance structure
  S_dyn_list    <- vector("list", T)                           # variance coming from crosslaggs, autoreg and innovations
  Sigma_e_list  <- vector("list", T)                           # innovations covariance
  S_target_list <- vector("list", T)                           # target covariance at observed level
  c_base_vec    <- numeric(T)                                  # the baseline covariance if residual covariance is 0 (meaning all coming from system + confounders)
  c_total_vec   <- numeric(T)                                  # the total covariance at observed level (including extra rho)

  # computing the variance structure at each wave
  # for each wave t from 1 to T
  for (t in 1:T) {

    # get the D matrix for this wave
    D_t <- D_list[[t]]

    # the confounder induced variance covariance at this wave is D_t Psi D_t'
    S_U_t <- D_t %*% Psi %*% t(D_t)

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
    Sigma_e_t <- S_dyn_t - A %*% S_dyn_t %*% t(A)

    # ensure symmetry
    Sigma_e_t <- (Sigma_e_t + t(Sigma_e_t))/2

    # store in container
    Sigma_e_list[[t]] <- Sigma_e_t
  }

  # prepare data frame to hold simulated data
  df <- matrix(NA, nrow = N, ncol = 2*T + k)

  # set column names
  colnames(df) <- c(paste0("x", 1:T),
                    paste0("y", 1:T),
                    paste0("c", 1:k))

  # add confounders to dataframe
  df[, (2*T + 1):(2*T + k)] <- U

  # simulate the first wave, which is different because there are no lagged values yet
  Ddyn <- mvtnorm::rmvnorm(
    n     = N,                                                   # number of individuals
    mean  = c(0, 0),                                             # mean 0 for X and Y
    sigma = S_dyn_list[[1]]                                      # variance covariance matrix at wave 1
  ) 

  # add the direct confounder effects
  obs1 <- Ddyn + U %*% t(D_list[[1]])

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
    obs <- Ddyn + U %*% t(D_list[[t]])

    # store
    df[, paste0("x", t)] <- obs[, 1]
    df[, paste0("y", t)] <- obs[, 2]
  }

  # return the simulated data frame
  return(df)
}
