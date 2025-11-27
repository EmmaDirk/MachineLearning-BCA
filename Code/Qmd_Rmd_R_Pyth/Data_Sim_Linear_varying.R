set.seed(427)

library(mvtnorm)
library(psych)
library(lavaan)
library(tidyverse)

N <- 100000
T <- 5

ax <- 0.25
ay <- 0.25
bx <- 0.10
by <- 0.10
rho <- 0.30

A <- matrix(c(ax, by, bx, ay), nrow=2, byrow=TRUE)

k <- 3
gamma_x <- c(0.30, 0.20, 0.10)
gamma_y <- c(0.30, 0.20, 0.10)
tau2 <- c(1,1,1)

B <- rbind(gamma_x, gamma_y)
Psi <- diag(tau2)

# --- time-varying effects (PERTURBATIONS AROUND BASELINE) ---
mean_change_x_vec <- c(0.00, 0.00, -0.30, 0.05, 0.05)
mean_change_y_vec <- c(0.00, 0.00, -0.20, 0.03, 0.03)
sd_change <- 0.005

B_list <- vector("list", T)
B_list[[1]] <- B   # baseline coefficients

for(t in 2:T){
  
  # mean changes for this wave (these shift the perturbation mean)
  mean_change_x <- mean_change_x_vec[t]
  mean_change_y <- mean_change_y_vec[t]
  
  # draw proportional perturbations
  change_x <- rnorm(k, mean = mean_change_x, sd = sd_change)
  change_y <- rnorm(k, mean = mean_change_y, sd = sd_change)
  
  # new coefficients = baseline * (1 + perturbation)
  new_gamma_x <- B_list[[1]][1, ] * (1 + change_x)
  new_gamma_y <- B_list[[1]][2, ] * (1 + change_y)
  
  B_list[[t]] <- rbind(new_gamma_x, new_gamma_y)
}
# ------------------------------------------------------------

find_c <- function(A, S_U, rho){

  f <- function(c){
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         2, 2, byrow=TRUE)
    S_dyn_c <- S_target_c - S_U
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    v1 <- Sigma_e_c[1,1]
    v2 <- Sigma_e_c[2,2]
    cov12 <- Sigma_e_c[1,2]
    corr_e <- cov12 / sqrt(v1*v2)
    corr_e - rho
  }
  
  uniroot(f, interval=c(-0.99, 0.99))$root
}

S_target_list <- vector("list", T)
S_U_list      <- vector("list", T)
S_dyn_list    <- vector("list", T)
Sigma_e_list  <- vector("list", T)

for(t in 1:T){

  B_t <- B_list[[t]]
  S_U_t <- B_t %*% Psi %*% t(B_t)
  S_U_list[[t]] <- S_U_t
  
  c_t <- find_c(A, S_U_t, rho)
  
  S_target_t <- matrix(c(1, c_t,
                         c_t, 1),
                       2, 2, byrow=TRUE)
  S_target_list[[t]] <- S_target_t
  
  S_dyn_t <- S_target_t - S_U_t
  S_dyn_t <- (S_dyn_t + t(S_dyn_t))/2
  S_dyn_list[[t]] <- S_dyn_t
  
  Sigma_e_t <- S_dyn_t - t(A) %*% S_dyn_t %*% A
  Sigma_e_t <- (Sigma_e_t + t(Sigma_e_t))/2
  Sigma_e_list[[t]] <- Sigma_e_t
}

Sigma_e1 <- S_dyn_list[[1]]
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

U <- rmvnorm(N, mean=rep(0,k), sigma=Psi)

Ddyn <- rmvnorm(N, mean=c(0,0), sigma=Sigma_e1)

df <- matrix(NA, nrow=N, ncol=2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  paste0("c", 1:k))

df[, (2*T+1):(2*T+k)] <- U

obs1 <- Ddyn + U %*% t(B_list[[1]])

df[,1]   <- obs1[,1]
df[,1+T] <- obs1[,2]

for(i in 2:T){
  
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma=Sigma_e_list[[i]])
  obs <- Ddyn + U %*% t(B_list[[i]])
  
  df[, i]   <- obs[,1]
  df[, i+T] <- obs[,2]
}

round(apply(df, 2, var), 3)
sigma <- cov(df)
round(sigma, 3)
describe(df)

df_sample <- df[sample(1:nrow(df), 10000), ]
lm_fit <- lm(x1 ~ c1, data = as.data.frame(df_sample))
residuals <- lm_fit$residuals

ggplot(df_sample, aes(x=c1, y=residuals))+
  geom_point(alpha=.5)+
  geom_smooth(method="loess", se=FALSE)+
  theme_minimal()
