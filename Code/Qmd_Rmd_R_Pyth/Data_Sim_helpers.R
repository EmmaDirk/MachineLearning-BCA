
# for starters we need a function to compute the total strength of the confounders
# now since we know that our confounders are unrelated, and their variances one
# we can compute the total effect size (or variance explained) as B*I*B^T = B*B^T = sum of squares of the elements of B

# B: 2 x k matrix of confounder effects
#     row 1: effects on X (gamma_x)
#     row 2: effects on Y (gamma_y)
#
# Returns:
# - S_U: 2x2 confounder-induced var-cov matrix
# - var_X: variance in X explained by confounders
# - var_Y: variance in Y explained by confounders
# - cov_XY: covariance between X and Y explained by confounders

confounder_strength <- function(B) {
  
  if (!is.matrix(B) || nrow(B) != 2) {
    stop("B must be a 2 x k matrix (rows = X, Y; columns = confounders).")
  }
  
  # with Psi = I, S_U = B %*% t(B)
  S_U <- B %*% t(B)
  
  var_X <- S_U[1, 1]
  var_Y <- S_U[2, 2]
  cov_XY <- S_U[1, 2]
  
  list(
    S_U   = S_U,
    var_X = var_X,
    var_Y = var_Y,
    cov_XY = cov_XY
  )
}

# now we want to create a function that computes a simple index to quantify the nonlinearity of the confounders
# which here is simply the nonlinear variance divided by total variance
# or put otherwise: B_nl*I*B_nl^T / (B*I*B^T) = B_nl*B_nl^T / (B*B^T)
# where B_nl is the matrix of nonlinear confounder effects

nonlinearity_index <- function(gamma_x, gamma_y, k_lin) {
  
  # helper
  eta <- function(g_lin, g_nl) {
    var_L <- sum(g_lin^2)
    var_N <- sum(g_nl^2)
    var_N / (var_L + var_N)
  }
  
  # split into linear and nonlinear parts
  k_total <- length(gamma_x)
  gamma_x_lin <- gamma_x[1:k_lin]
  gamma_y_lin <- gamma_y[1:k_lin]
  
  gamma_x_nl  <- gamma_x[(k_lin + 1):k_total]
  gamma_y_nl  <- gamma_y[(k_lin + 1):k_total]
  
  # compute nonlinearity for X and Y
  eta_X <- eta(gamma_x_lin, gamma_x_nl)
  eta_Y <- eta(gamma_y_lin, gamma_y_nl)
  
  # return their average
  mean(c(eta_X, eta_Y))
}

# test it
nonlinearity_index(gamma_x, gamma_y, k_lin)

# lets now produce a few different gamma matrixes and test the function
# mostly linear
gamma_x1 <- c(0.30, 0.20, 0.10,
               0.05, 0.02, 0.00,
               0.02, 0.01, 0.00)

gamma_y1 <- c(0.30, 0.20, 0.10,
               0.05, 0.02, 0.00,
               0.02, 0.01, 0.00)

nonlinearity_index(gamma_x1, gamma_y1, k_lin)  # expect low value

# now note that this could be expanded upon by simply adding the covariance matrix
# psi of the confounders if we want our confounders to be related.   
# below we work out the generalised versions