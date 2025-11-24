# the purpose of this script is to understand what happens if you tryt to sample from a thin surface,
# i.e., a space where the one dimension is scalar.


### 2D circle inside square ###
set.seed(1)

# number of random points in the square
n <- 5000

# sample uniformly from [-1,1]^2
theta1 <- runif(n, -1, 1)
theta2 <- runif(n, -1, 1)

# radius of the circle (constraint)
r <- 0.8

# distance from origin
dist2 <- sqrt(theta1^2 + theta2^2)

# points "near" the circle: thin band of width eps
eps <- 0.02
on_circle <- abs(dist2 - r) <= eps

# plot
plot(theta1, theta2,
     col = ifelse(on_circle, "red", "grey80"),
     pch = 16, cex = 0.6,
     xlab = expression(theta[1]),
     ylab = expression(theta[2]),
     main = "Thin surface in 2D: circle inside square")
abline(h = 0, v = 0, col = "grey50", lty = 3)
box()


### 3D sphere inside cube ###
library(scatterplot3d)

set.seed(2)

n <- 8000

# sample uniformly from [-1,1]^3
theta1 <- runif(n, -1, 1)
theta2 <- runif(n, -1, 1)
theta3 <- runif(n, -1, 1)

# sphere radius
r <- 0.9

# distance from origin
dist3 <- sqrt(theta1^2 + theta2^2 + theta3^2)

# narrow band around the sphere
eps <- 0.03
on_sphere <- abs(dist3 - r) <= eps

# 3D scatter: cube in grey, thin shell in red
scatterplot3d(theta1, theta2, theta3,
              color = ifelse(on_sphere, "red", "grey80"),
              pch   = 16,
              xlab  = expression(theta[1]),
              ylab  = expression(theta[2]),
              zlab  = expression(theta[3]),
              main  = "Thin surface in 3D: sphere shell inside cube")

# now if we want to sample from such a space, the probability of being "on the surface" is very low
# a solution to this is to make the scalar dimension have a little more thickness by accepting some
# tolerance around the constraint. It would look as follows:

### 2D circle with thickness ###

set.seed(1)

# number of random points in the square
n <- 5000

# sample uniformly from [-1,1]^2
theta1 <- runif(n, -1, 1)
theta2 <- runif(n, -1, 1)

# radius of the circle (constraint)
r <- 0.8

# distance from origin
dist2 <- sqrt(theta1^2 + theta2^2)

# different tolerances to try
eps_values <- c(0.005, 0.02, 0.08)

par(mfrow = c(1, length(eps_values)))  # 1 row, 3 columns

for (eps in eps_values) {
  
  on_circle <- abs(dist2 - r) <= eps
  
  plot(theta1, theta2,
       col = ifelse(on_circle, "red", "grey80"),
       pch = 16, cex = 0.6,
       xlab = expression(theta[1]),
       ylab = expression(theta[2]),
       main = paste0("eps = ", eps))
  
  # draw the circle itself (for reference)
  ang <- seq(0, 2*pi, length.out = 200)
  lines(r * cos(ang), r * sin(ang), col = "blue", lwd = 2)
  
  box()
}
par(mfrow = c(1,1))

### 3D sphere with thickness ###

set.seed(2)

n <- 8000

# sample uniformly from [-1,1]^3
theta1 <- runif(n, -1, 1)
theta2 <- runif(n, -1, 1)
theta3 <- runif(n, -1, 1)

# sphere radius
r <- 0.9

# distance from origin
dist3 <- sqrt(theta1^2 + theta2^2 + theta3^2)

# different tolerances
eps_values <- c(0.01, 0.04, 0.12)

par(mfrow = c(1, length(eps_values)))

for (eps in eps_values) {
  
  on_sphere <- abs(dist3 - r) <= eps
  
  scatterplot3d(theta1, theta2, theta3,
                color = ifelse(on_sphere, "red", "grey80"),
                pch   = 16,
                xlab  = expression(theta[1]),
                ylab  = expression(theta[2]),
                zlab  = expression(theta[3]),
                main  = paste0("3D shell, eps = ", eps))
}

par(mfrow = c(1,1))
