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

# nonlinear (sinusoidal) transformation parameters
omega <- 3

# simulate base confounders: C1..C3 ~ N(0,1), independent
Psi_base <- diag(rep(1, k_lin))                            
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    
                 sigma = Psi_base)

# build sinusoidal transforms S_j = scaled sin(omega * C_j)
U_sin <- apply(U_lin, 2, function(col) {

  z_raw <- sin(omega * col)

  z_perp <- resid(lm(z_raw ~ col))

  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

colnames(U_lin) <- paste0("c", 1:k_lin)                    # C1, C2, C3
colnames(U_sin) <- paste0("c", 1:k_lin, "_sin")            # S1, S2, S3

# combine linear and sinusoidal confounders
U <- cbind(U_lin, U_sin)                                   # N x k, where k = 2 * k_lin
k <- ncol(U)

Psi <- cov(U)

# effects for the sinusoidal parts on X and Y
gamma_x_sin <- c(0.15, 0.10, 0.05)
gamma_y_sin <- c(0.15, 0.10, 0.05)

# order of entries in gamma_x / gamma_y corresponds to (C1, C2, C3, S1, S2, S3)
gamma_x <- c(gamma_x_lin, gamma_x_sin)
gamma_y <- c(gamma_y_lin, gamma_y_sin)

### --- NEW: handle B or B_list --- ###

# baseline 2 x k effect matrix from all confounders to (X, Y)
B_single <- rbind(gamma_x, gamma_y)

# If the user has defined a time-varying B_list (list of length T with 2 x k matrices),
# we use that. Otherwise we fall back to a constant B over time.
if (!exists("B_list")) {
  # fallback: same B at every wave
  B_list <- replicate(T, B_single, simplify = FALSE)
}

# compute average confounder-induced variance across waves:
# S_U_t = B_t Psi B_t'
S_U_list <- lapply(B_list, function(Bt) Bt %*% Psi %*% t(Bt))
S_U <- Reduce("+", S_U_list) / length(S_U_list)

### varcov matrices ###
# we want the variance of each variable to be 1
S_target <- diag(2)                                       

# function to find the required covariance between X and Y
# such that the innovation correlation is rho

find_c <- function(A, S_U, rho) {
  f <- function(c) {

    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    
    S_dyn_c <- S_target_c - S_U
    
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    
    v1 <- Sigma_e_c[1, 1]
    v2 <- Sigma_e_c[2, 2]
    cov12 <- Sigma_e_c[1, 2]
    corr_e <- cov12 / sqrt(v1 * v2)
    
    corr_e - rho
  }

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
eS[eS < 1e-10]

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

# --- NEW: use B_list[[1]] for wave 1 ---
obs1 <- Ddyn + U %*% t(B_list[[1]])

df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# subsequent waves
for(i in 2:T){ 
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 
  
  # --- NEW: use B_list[[i]] for wave i ---
  obs  <- Ddyn + U %*% t(B_list[[i]]) 
  
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2] 
}

df <- as.data.frame(df)

# checks
round(apply(df, 2, var), 3)
sigma <- cov(df)
round(sigma, 3)
describe(df)

# quick residual check for nonlinearity
df_sample <- df[sample(1:nrow(df), 10000), ]
lm_fit <- lm(x1 ~ c1, data = df_sample)
residuals <- lm_fit$residuals

ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# write out the data
# write.csv(df, file = "Thesis/Code/Data/Data_Sim_Sin_Pop.csv", row.names = FALSE)
