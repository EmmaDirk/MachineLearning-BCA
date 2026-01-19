## ============================================================
##  CLPM with linear/quadratic/cubic confounders
##  AND time-varying effects of confounders (B_t)
##  while keeping stationarity and Var(X_t)=Var(Y_t)=1
## ============================================================

set.seed(427)

library(mvtnorm)
library(psych)
library(lavaan)
library(tidyverse)

### CLPM parameters ###
N <- 100000          # number of people
T <- 5               # number of time points

ax  <- 0.25          # X_(t-1) -> X_t
ay  <- 0.25          # Y_(t-1) -> Y_t
bx  <- 0.10          # Y_(t-1) -> X_t
by  <- 0.10          # X_(t-1) -> Y_t
rho <- 0.30          # residual correlation of innovations

# path matrix A
A <- matrix(c(ax, by,
              bx, ay),
            nrow = 2, byrow = TRUE)

### confounder parameters ###
k_lin   <- 3         # number of base/linear confounders
k_quad  <- 3         # number of quadratic confounders
k_cubic <- 3         # number of cubic confounders
k       <- k_lin + k_quad + k_cubic   # total number of confounders

# base effect sizes (for t = 1, or "baseline" level)
gamma_x <- c(0.30, 0.20, 0.10,   # linear
             0.05, 0.02, 0.00,   # quadratic
             0.02, 0.10, 0.00)   # cubic

gamma_y <- c(0.30, 0.20, 0.10,   # linear
             0.05, 0.02, 0.00,   # quadratic
             0.02, 0.10, 0.00)   # cubic

## ------------------------------------------------------------
## Time-varying B matrices (effects of confounders vary by wave)
## ------------------------------------------------------------

# scaling factors for each wave (you can change this pattern)
scales <- c(1.00, 0.90, 0.80, 0.70, 0.60)

# list of B_t (2 x k each), one for each time point
B_list <- lapply(scales, function(s) {
  rbind(gamma_x * s,
        gamma_y * s)
})

### generate confounders (time-invariant) ###

# variances of base confounders: N(0,1), independent
tau2      <- rep(1, k_lin)
Psi_base  <- diag(tau2)

# linear confounders
U_lin <- rmvnorm(N, mean = rep(0, k_lin), sigma = Psi_base)

# quadratic confounders: orthogonalized, Var = 1
U_quad <- apply(U_lin, 2, function(col) {
  z_raw  <- col^2
  z_perp <- resid(lm(z_raw ~ col))
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

# cubic confounders: orthogonalized, Var = 1
U_cubic <- apply(U_lin, 2, function(col) {
  z_raw  <- col^3
  z_perp <- resid(lm(z_raw ~ col))
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

# full confounder matrix: N x k
U <- cbind(U_lin, U_quad, U_cubic)

# empirical covariance matrix of all confounder terms
Psi <- cov(U)

### stationarity calculations with time-varying B_t ###

# target Var(X_t) = Var(Y_t) = 1, unknown Cov(X_t, Y_t) = c
S_target <- diag(2)

# confounder-induced variance per wave: S_U_t = B_t * Psi * B_t'
S_U_list <- lapply(B_list, function(B_t) B_t %*% Psi %*% t(B_t))

# use the average confounder contribution over waves to enforce stationarity
S_U <- Reduce("+", S_U_list) / length(S_U_list)

# function to find stationary covariance c such that innovation corr = rho
find_c <- function(A, S_U, rho) {
  f <- function(c) {
    # stationary target covariance of (X_t, Y_t)
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    # dynamic part (exclude confounder variance)
    S_dyn_c <- S_target_c - S_U

    # innovation covariance required for stationarity
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A

    v1    <- Sigma_e_c[1, 1]
    v2    <- Sigma_e_c[2, 2]
    cov12 <- Sigma_e_c[1, 2]
    corr_e <- cov12 / sqrt(v1 * v2)

    corr_e - rho
  }

  uniroot(f, interval = c(-0.99, 0.99))$root
}

# find stationary covariance c between X and Y
c_stat <- find_c(A, S_U, rho)

# overwrite S_target with this covariance
S_target <- matrix(c(1, c_stat,
                     c_stat, 1),
                   nrow = 2, byrow = TRUE)

# dynamic part variance (what's left after confounders)
S_dyn <- S_target - S_U

# check positive semidefiniteness
eS <- eigen((S_dyn + t(S_dyn)) / 2)$values
if (any(eS < 1e-10)) {
  warning("S_dyn has (near-)negative eigenvalues; check parameters.")
}

# innovation covariance for t >= 2
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e)) / 2

# innovation covariance for t = 1 (no lag yet)
# S_target = S_U + S_dyn  =>  Sigma_e1 = S_dyn
Sigma_e1 <- (S_dyn + t(S_dyn)) / 2

### simulate data ###

# dynamic part for t = 1
Ddyn <- rmvnorm(N,
                mean  = c(0, 0),
                sigma = Sigma_e1)

# empty matrix: 2*T for x1..xT, y1..yT, plus k confounders
df <- matrix(NA, nrow = N, ncol = 2 * T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  paste0("c", 1:k))

# store confounders
df[, (2 * T + 1):(2 * T + k)] <- U

# first wave: add confounder effects B_1
obs1 <- Ddyn + U %*% t(B_list[[1]])

df[, 1]    <- obs1[, 1]     # x1
df[, 1 + T] <- obs1[, 2]    # y1

# subsequent waves t = 2..T
for (i in 2:T) {
  # dynamic part
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e)

  # add time-specific confounder effects B_i
  obs <- Ddyn + U %*% t(B_list[[i]])

  df[, i]     <- obs[, 1]      # x_t
  df[, i + T] <- obs[, 2]      # y_t
}

df <- as.data.frame(df)

## quick checks
round(apply(df, 2, var), 3)      # variances ≈ 1
sigma <- cov(df)
round(sigma[1:(2*T), 1:(2*T)], 3)  # var-cov of X/Y only

describe(df)

# optional: save
# write.csv(df, file = "Thesis/Code/Data/Data_Sim_Cubic_TimeVaryingB_Pop.csv",
#           row.names = FALSE)
