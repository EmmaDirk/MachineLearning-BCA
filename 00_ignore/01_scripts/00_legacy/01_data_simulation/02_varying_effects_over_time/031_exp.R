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

# curvature parameter for the exponential (controls how fast it explodes)
b_exp <- 0.7                                               # you can tune this

# simulate base confounders: C1..C3 ~ N(0,1), independent
Psi_base <- diag(rep(1, k_lin))                            
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    
                 sigma = Psi_base)                         

colnames(U_lin) <- paste0("c", 1:k_lin)                    # C1, C2, C3

# exponential-transformed nonlinear confounders:
# E_j = standardized exp(b_exp * C_j), orthogonalized to C_j
U_exp <- apply(U_lin, 2, function(col) {
  z_raw <- exp(b_exp * col)                                # skewed, heavy right tail
  
  # remove the linear part of col (orthogonalize wrt the linear confounder)
  z_perp <- resid(lm(z_raw ~ col))
  
  # standardize to mean 0 and var 1
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

colnames(U_exp) <- paste0("c", 1:k_lin, "_exp")            # E1, E2, E3

# combine linear and exponential confounders
U <- cbind(U_lin, U_exp)                                   # N x k, where k = 2 * k_lin
k <- ncol(U)                                               # total number of confounder variables

# empirical varcov matrix of (C, exp(C)) confounders
Psi <- cov(U)

# effects for the exponential parts on X and Y (usually smaller)
gamma_x_exp <- c(0.10, 0.07, 0.05)
gamma_y_exp <- c(0.10, 0.07, 0.05)

# order of entries corresponds to (C1, C2, C3, E1, E2, E3)
gamma_x <- c(gamma_x_lin, gamma_x_exp)
gamma_y <- c(gamma_y_lin, gamma_y_exp)

# baseline 2 x k effect matrix
B_single <- rbind(gamma_x, gamma_y)

### allow time-varying B via B_list ###

# If B_list exists (list length T with 2 x k matrices), use that.
# Otherwise, default to constant B across time.
if (!exists("B_list")) {
  B_list <- replicate(T, B_single, simplify = FALSE)
}

### varcov matrices ###

# we want the variance of each variable to be 1
S_target <- diag(2)                                        # initial 2x2 identity (var = 1, cov = 0)

# variance contributed by confounders over time:
# S_U_t = B_t * Psi * B_t'
S_U_list <- lapply(B_list, function(Bt) Bt %*% Psi %*% t(Bt))
S_U  <- Reduce("+", S_U_list) / length(S_U_list)           # mean across time

# find needed covariance c between X and Y 
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
    
    # compute correlation of innovations
    v1 <- Sigma_e_c[1, 1]
    v2 <- Sigma_e_c[2, 2]
    cov12 <- Sigma_e_c[1, 2]
    corr_e <- cov12 / sqrt(v1 * v2)
    
    corr_e - rho
  }

  # find c such that f(c) = 0
  uniroot(f, interval = c(-0.99, 0.99))$root
}

# required covariance between X and Y
c_stat <- find_c(A, S_U, rho)

# stationary covariance of (X, Y)
S_target <- matrix(c(1, c_stat,
                     c_stat, 1),
                   nrow = 2, byrow = TRUE)

# variance by the dynamic process:
S_dyn <- S_target - S_U
S_dyn <- (S_dyn + t(S_dyn))/2                             # symmetrise

# check positive semidefinite
eS <- eigen(S_dyn)$values
eS[eS < 1e-10]

# innovation covariance matrix
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2

# For t = 1, target: S_target = S_U + Sigma_e1
Sigma_e1 <- S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

# t = 1: no lag yet, so Ddyn = e1 with Sigma_e1
Ddyn <- rmvnorm(N,
  mean  = c(0,0),
  sigma = Sigma_e1)

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  colnames(U))

# save the confounders (C1..C3 and E1..E3)
df[, (2*T + 1):(2*T + k)] <- U

# first wave: use B_list[[1]]
obs1 <- Ddyn + U %*% t(B_list[[1]])

df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# remaining waves: time-varying B_list[[i]]
for(i in 2:T){ 
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 
  obs  <- Ddyn + U %*% t(B_list[[i]]) 
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2] 
}

df <- as.data.frame(df)

# basic checks
round(apply(df, 2, var), 3)
sigma <- cov(df)
round(sigma, 3)
describe(df)

# quick residual check for nonlinearity (just one example)
df_sample <- df[sample(1:nrow(df), 10000), ]
lm_fit <- lm(x1 ~ c1, data = df_sample)
residuals <- lm_fit$residuals

ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# write.csv(df, file = "Thesis/Code/Data/Data_Sim_exponential_Blist_Pop.csv", row.names = FALSE)
