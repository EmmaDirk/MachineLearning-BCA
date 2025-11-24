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

# curvature parameter: lager = quicker "turning" around a point
# meaning a stronger non-linear effect.
a_tanh <- 2

# simulate base confounders: C1..C3 ~ N(0,1), independent
Psi_base <- diag(rep(1, k_lin))                            # variance-covariance matrix of base confounders (I in this case)
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    # mean of confounders = 0, since we want standardised effect sizes
                 sigma = Psi_base)                         # C1..C3 ~ N(0,1), independent

# name the confounder columns
colnames(U_lin) <- paste0("c", 1:k_lin)                    # C1, C2, C3

# create tanh-transformed nonlinear confounders:
# T_j = standardized tanh(a_tanh * C_j), so each has mean ~0 and var ~1
# we also remove the linear part of C_j from T_j, so no dependency remains between the two (orthogonalization).
# this will later be useful to separate linear and non-linear effects
# when we want to calculate a measure of non-linearity later.
U_tanh <- apply(U_lin, 2, function(col) {

  # raw tanh transform of each confounder
  z_raw <- tanh(a_tanh * col)

  # remove the linear part of col (orthogonalize wrt the linear confounder)
  z_perp <- resid(lm(z_raw ~ col))

  # standardize to mean 0 and var 1
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

# name the confounder columns
colnames(U_tanh) <- paste0("c", 1:k_lin, "_tanh")          # T1, T2, T3

# combine linear and tanh confounders into one matrix U
# columns: first the k_lin linear confounders, then k_lin tanh confounders
U <- cbind(U_lin, U_tanh)                                  # N x k, where k = 2 * k_lin
k <- ncol(U)                                               # total number of confounder variables

# empirical varcov matrix of (C, tanh(C)) confounders
Psi <- cov(U)

# effects for the tanh parts on X and Y
# you can tune these (often smaller than linear effects)
gamma_x_tanh <- c(0.15, 0.10, 0.05)
gamma_y_tanh <- c(0.15, 0.10, 0.05)

# order of entries in gamma_x / gamma_y corresponds to the column order in U:
# (C1, C2, C3, T1, T2, T3)
gamma_x <- c(gamma_x_lin, gamma_x_tanh)
gamma_y <- c(gamma_y_lin, gamma_y_tanh)

# we can now build matrix B, which is 2 x k
B <- rbind(gamma_x, gamma_y)

### varcov matrixes ###
# we want the variance of each variable to be 1
S_target <- diag(2)                                        # initial 2x2 identity (var = 1, cov = 0)

# variance contributed by confounders
# S_U = B * Psi * B'
S_U  <- B %*% Psi %*% t(B)

# find covariance c between X and Y needed for target residual correlation rho
find_c <- function(A, S_U, rho) {

  f <- function(c) {
    # target stationary covariance of (X_t, Y_t)
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    
    # dynamic part variance (remove confounder-induced variance)
    S_dyn_c <- S_target_c - S_U
    
    # innovation covariance matrix under stationarity
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    
    # correlation of innovations
    v1 <- Sigma_e_c[1, 1]
    v2 <- Sigma_e_c[2, 2]
    cov12 <- Sigma_e_c[1, 2]
    corr_e <- cov12 / sqrt(v1 * v2)
    
    corr_e - rho
  }

  uniroot(f, interval = c(-0.99, 0.99))$root
}

# required covariance between X and Y
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

# For t = 1, there is no carryover from a previous wave. The observed var-cov should be:
# Var([X1,Y1]) = S_U + Sigma_e1. To keep the same S_target at t=1, set:
Sigma_e1 <- S_dyn                             # because S_target = S_U + S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

# simulate panel data 

#  for t = 1: no lag yet, so Ddyn = e1 with Sigma_e1
Ddyn <- rmvnorm(N,                            
  mean  = c(0,0),                              
  sigma = Sigma_e1)                            

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  colnames(U))

# save the confounders (C1..C3 and T1..T3)
df[, (2*T + 1):(2*T + k)] <- U

# first wave: add direct confounder effects
obs1 <- Ddyn + U %*% t(B)

df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# remaining waves
for(i in 2:T){ 
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 
  obs  <- Ddyn + U %*% t(B) 
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2] 
}

# basic checks
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

write.csv(df, file = "Thesis/Code/Data/Data_Sim_tanh_Pop.csv", row.names = FALSE)
