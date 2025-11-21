# reproducibility
set.seed(427) 

# load required libraries
library(mvtnorm)
library(psych)                                             # for descriptives
library(lavaan)                                            # for fitting CLPM
library(tidyverse)

### CLPM parameters ###
N <- 100000                                                # number of people
T <- 5                                                     # number of time points

ax <- 0.25                                                 # effect of X_(t-1) on X_t
ay <- 0.25                                                 # effect of Y_(t-1) on Y_t
bx <- 0.10                                                 # effect of Y_(t-1) on X_t
by <- 0.10                                                 # effect of X_(t-1) on Y_t
rho <- 0.30                                                # residual correlation of innovations
                                                          
# create the path matrix A from these effects
A <- matrix(c(ax, by, 
              bx, ay), 
            nrow = 2, byrow = TRUE)

### confounder parameters ###
k_lin <- 3                                                 # number of "linear" confounders (baseline C1..C3)

# baseline linear effects of confounders on X_t and Y_t
gamma_x_lin <- c(0.30, 0.20, 0.10)                         # effect of confounders on X_t
gamma_y_lin <- c(0.30, 0.20, 0.10)                         # effect of confounders on Y_t

# we are now going to add a nonlinear (sinusoidal) transformation of the confounders
# to simulate a non-linear confounding effect
# frequency of the oscillation (how fast the sine wave oscillates)
omega <- 3

# simulate base confounders: C1..C3 ~ N(0,1), independent
Psi_base <- diag(rep(1, k_lin))                            # variance-covariance matrix of base confounders (I in this case)
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    # mean of confounders = 0, since we want standardised effect sizes
                 sigma = Psi_base)

# build sinusoidal transforms S_j = scaled sin(omega * C_j)
# so each has mean ~0 and var ~1
# we also remove the linear part of C_j from S_j, so no dependency remains between the two (orthogonalization).
# this will later be useful to separate linear and non-linear effects
# when we want to calculate a measure of non-linearity later.
U_sin <- apply(U_lin, 2, function(col) {

  # raw sinusoidal transform of each confounder
  z_raw <- sin(omega * col)

  # remove the linear part of col (orthogonalize wrt the linear confounder)
  z_perp <- resid(lm(z_raw ~ col))

  # scale the variables such that they have mean 0 and variance 1
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})


# name the confounder columns to distinguish linear and sinusoidal parts
colnames(U_lin) <- paste0("c", 1:k_lin)                    # C1, C2, C3
colnames(U_sin) <- paste0("c", 1:k_lin, "_sin")            # S1, S2, S3 (nonlinear parts)

# combine linear and sinusoidal confounders into one matrix U
# columns: first the k_lin linear confounders, then k_lin sinusoidal confounders
U <- cbind(U_lin, U_sin)                                   # N x k, where k = 2 * k_lin

# total number of confounder variables (C and sinusoid versions)
k <- ncol(U)

# empirical var-cov matrix of all confounders
Psi <- cov(U)

# effects for the sinusoidal parts on X and Y
# larger effects mean: greater non-linearity in the confounding
gamma_x_sin <- c(0.15, 0.10, 0.05)
gamma_y_sin <- c(0.15, 0.10, 0.05)

# order of entries in gamma_x / gamma_y corresponds to the column order in U:
# (C1, C2, C3, S1, S2, S3)
gamma_x <- c(gamma_x_lin, gamma_x_sin)
gamma_y <- c(gamma_y_lin, gamma_y_sin)

# 2 x k effect matrix from all confounders to (X, Y)
B <- rbind(gamma_x, gamma_y)

### varcov matrixes ###
# we want the variance of each variable to be 1
S_target <- diag(2)                                        # since we have two variables, X and Y, we create a 2x2 identity matrix. This assumes NO covariance between X and Y
                                                           # further than the effects that were specified before

# variance contributed by confounders:
# S_U = B * Psi * B'
S_U  <- B %*% Psi %*% t(B)

# function to find the required covariance between X and Y
# such that the innovation correlation is rho

find_c <- function(A, S_U, rho) {

  # function of c: difference between actual innovation correlation and target rho
  # we need to calculate how much covariance is induced by c and the lagged relationships
  f <- function(c) {

    # target stationary covariance of (X_t, Y_t)
    # the variances are 1 since we specified them to be
    # and the covariance is c, which we are searching for
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    
    # dynamic part variance
    # first we get rid of the confounder-induced variance
    S_dyn_c <- S_target_c - S_U
    
    # now we can compute the innovation covariance matrix
    # based on the stationarity constraints
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    
    # compute correlation of innovations
    # the variance of innovations of X is v1
    v1 <- Sigma_e_c[1, 1]
    # the variance of innovations of Y is v2
    v2 <- Sigma_e_c[2, 2]
    # their covariance is cov12
    cov12 <- Sigma_e_c[1, 2]
    # their correlation is then:
    corr_e <- cov12 / sqrt(v1 * v2)
    
    # this is then the difference between the desired correlation and the actual correlation
    corr_e - rho
  }

  # now this function f(c) should be zero when we have found the correct c
  # as such we can look for the value of c that makes f(c) = 0.
  uniroot(f, interval = c(-0.99, 0.99))$root
}

# use the function to find the required covariance between X and Y
c_stat <- find_c(A, S_U, rho)

# overwrite S_target to include the stationary covariance between X and Y
S_target <- matrix(c(1, c_stat,
                     c_stat, 1),
                   nrow = 2, byrow = TRUE)

# variance by the dynamic process:
S_dyn <- S_target - S_U

# check positive semidefinite
eS <- eigen((S_dyn + t(S_dyn))/2)$values
eS[eS<1e-10]

# innovation covariance matrix
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2

# For t = 1, no carryover
Sigma_e1 <- S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

#  for t = 1: no lag yet, so Ddyn = e1 with Sigma_e1
Ddyn <- rmvnorm(N,                            
  mean  = c(0,0),                              
  sigma = Sigma_e1)                            

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  colnames(U))

# save the confounders (C1..C3 and S1..S3)
df[, (2*T + 1):(2*T + k)] <- U

# first wave: add direct confounder effects
obs1 <- Ddyn + U %*% t(B)

df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# subsequent waves
for(i in 2:T){ 
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 
  obs  <- Ddyn + U %*% t(B) 
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2] 
}

# checks
round(apply(df, 2, var), 3)
sigma <- cov(df)
round(sigma, 3)
describe(df)

# quick residual check for nonlinearity
df_sample <- df[sample(1:nrow(df), 10000), ]
lm_fit <- lm(x1 ~ c1, data = as.data.frame(df_sample))
residuals <- lm_fit$residuals

ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()
