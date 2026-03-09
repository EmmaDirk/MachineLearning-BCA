# this script serves to test the delta sampler
#
#
# Recall that the sampler:
# - draws delta coefficients for X and Y at time t = 1,
# - enforces that the variance explained by the linear part equals Vlin_star,
# - enforces that the variance explained by the interaction part equals VNL_star,
# - enforces that the total variance explained equals R2_total,
# - optionally includes all 2-way and/or 3-way interactions,
# - enforces absolute bounds on each coefficient,
# - errors when logical constraints are violated.
#
# We therefore test:
#
# 1) That the sampler errors under logically inconsistent inputs:
#    a) R2_total outside [0, 1].
#    b) R2_interaction outside [0, 1].
#    c) R2_interaction > 0 but no interactions toggled on.
#    d) Interactions toggled on but R2_interaction = 0.
#    e) k too small for requested interactions.
#    f) Sigma missing or invalid.
#
# 2) That the returned object has the correct structure:
#    a) A list.
#    b) A numeric Delta matrix.
#    c) Two rows named "X" and "Y".
#    d) The correct number of columns given k and the toggles.
#    e) Correct column naming and ordering.
#    f) Omega and Omega_blocks are returned.
#
# 3) That the mathematical constraints are satisfied:
#    a) The variance explained by the linear part equals Vlin_star.
#    b) The variance explained by the interaction part equals VNL_star.
#    c) The total variance explained equals R2_total.
#
# 4) That min_abs / max_abs bounds are respected.
#
# 5) That infeasible constraints lead to an error.
#
# 6) That the sampler is reproducible under set.seed().
# ----------------------------------------------------------------------------------------------------------------

# libraries
library(here)

# source script
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))

# reproducibility
set.seed(1)

# ------------------------------------------------------------------------------------------------------------------
# helpers used in the tests
# ------------------------------------------------------------------------------------------------------------------

# example covariance matrix of the standardized base confounders
Sigma <- matrix(c(
  1.0, 0.3, 0.2,
  0.3, 1.0, 0.4,
  0.2, 0.4, 1.0
), nrow = 3, byrow = TRUE)

# helper to compute the implied variance decomposition for one row of Delta
# this is to test that the mathematical constraints are satisfied
get_variance_components <- function(delta_row, out, k, include_2way, include_3way) {

  # number of possible interaction terms
  m2 <- if (include_2way) choose(k, 2) else 0
  m3 <- if (include_3way) choose(k, 3) else 0

  # split the coefficient vector into the linear, 2way, and 3way parts
  dL <- delta_row[1:k]
  d2 <- if (m2 > 0) delta_row[(k + 1):(k + m2)] else numeric(0)
  d3 <- if (m3 > 0) delta_row[(k + m2 + 1):(k + m2 + m3)] else numeric(0)

  # extract Omega blocks
  Sigma_block   <- out$Omega_blocks$Sigma
  Omega22_block <- out$Omega_blocks$Omega22
  Omega33_block <- out$Omega_blocks$Omega33
  OmegaL3_block <- out$Omega_blocks$OmegaL3

  # linear variance contribution
  A <- as.numeric(t(dL) %*% Sigma_block %*% dL)

  # 2way interaction variance contribution
  B <- if (m2 > 0) as.numeric(t(d2) %*% Omega22_block %*% d2) else 0

  # 3way interaction variance contribution
  C <- if (m3 > 0) as.numeric(t(d3) %*% Omega33_block %*% d3) else 0

  # covariance contribution between linear and 3way parts
  D <- if (m3 > 0) as.numeric(t(dL) %*% OmegaL3_block %*% d3) else 0

  # return the decomposition exactly as used by the sampler
  list(
    V_lin = A + D,
    V_int = B + C + D,
    V_tot = A + B + C + 2 * D,
    A = A,
    B = B,
    C = C,
    D = D
  )
}

# ------------------------------------------------------------------------------------------------------------------
# 1: test that the sampler errors under logically inconsistent inputs
# ------------------------------------------------------------------------------------------------------------------

# a) R2_total outside [0, 1]
# SHOULD ERROR
sample_delta_1(k = 3,
               Sigma = Sigma,
               R2_total = 1.1,
               R2_interaction = 0,
               include_2way = FALSE,
               include_3way = FALSE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# b) R2_interaction outside [0, 1]
# SHOULD ERROR
sample_delta_1(k = 3,
               Sigma = Sigma,
               R2_total = 0.5,
               R2_interaction = 1.1,
               include_2way = FALSE,
               include_3way = FALSE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# c) R2_interaction > 0 but no interactions toggled on
# SHOULD ERROR
sample_delta_1(k = 3,
               Sigma = Sigma,
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = FALSE,
               include_3way = FALSE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# d) Interactions toggled on but R2_interaction = 0
# SHOULD ERROR
sample_delta_1(k = 3,
               Sigma = Sigma,
               R2_total = 0.5,
               R2_interaction = 0,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# e) k too small for requested interactions
# SHOULD ERROR
sample_delta_1(k = 2,
               Sigma = diag(2),
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# f1) Sigma missing
# SHOULD ERROR
sample_delta_1(k = 3,
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# f2) Sigma not symmetric
# SHOULD ERROR
Sigma_nonsym <- matrix(c(
  1.0, 0.3, 0.2,
  0.1, 1.0, 0.4,
  0.2, 0.4, 1.0
), nrow = 3, byrow = TRUE)

sample_delta_1(k = 3,
               Sigma = Sigma_nonsym,
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# f3) Sigma diagonal not equal to 1
# SHOULD ERROR
Sigma_bad_diag <- matrix(c(
  2.0, 0.3, 0.2,
  0.3, 1.0, 0.4,
  0.2, 0.4, 1.0
), nrow = 3, byrow = TRUE)

sample_delta_1(k = 3,
               Sigma = Sigma_bad_diag,
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

# ------------------------------------------------------------------------------------------------------------------
# 2: test that the returned object has the correct structure
# ------------------------------------------------------------------------------------------------------------------

# first visually inspect the output
# should be:
# a) a list.
# b) a numeric Delta matrix with 2 rows and k + choose(k, 2) + choose(k, 3) columns.
# c) two row names "X" and "Y"
# d) the correct number of columns given k and the toggles.
# e) correct column naming and ordering.
# f) Omega and Omega_blocks are returned.
d <- sample_delta_1(k = 3,
                    Sigma = Sigma,
                    R2_total = 0.5,
                    R2_interaction = 0.5,
                    include_2way = TRUE,
                    include_3way = TRUE,
                    min_abs = 0,
                    max_abs = 1,
                    max_tries = 1000000)

# visually inspect
d

# check structure
str(d)

# inspect the Delta matrix directly
d$Delta

# inspect the full covariance matrix directly
d$Omega

# inspect the Omega blocks directly
d$Omega_blocks

# basic checks
is.list(d)
is.matrix(d$Delta)
is.numeric(d$Delta)
identical(rownames(d$Delta), c("X", "Y"))
ncol(d$Delta) == 7
identical(colnames(d$Delta), c("c1", "c2", "c3", "c1:2", "c1:3", "c2:3", "c1:2:3"))

# ------------------------------------------------------------------------------------------------------------------
# 3: test that the mathematical constraints are satisfied
# ------------------------------------------------------------------------------------------------------------------

# target variance decomposition
Vlin_star <- (1 - 0.5) * 0.5
VNL_star  <- 0.5 * 0.5

# compute the implied variance decomposition for X and Y
vx <- get_variance_components(delta_row = d$Delta["X", ],
                              out = d,
                              k = 3,
                              include_2way = TRUE,
                              include_3way = TRUE)

vy <- get_variance_components(delta_row = d$Delta["Y", ],
                              out = d,
                              k = 3,
                              include_2way = TRUE,
                              include_3way = TRUE)

# a) variance explained by the linear part equals Vlin_star
# SHOULD BE TRUE
abs(vx$V_lin - Vlin_star) < 1e-10
abs(vy$V_lin - Vlin_star) < 1e-10

# b) variance explained by the interaction part equals VNL_star
# SHOULD BE TRUE
abs(vx$V_int - VNL_star) < 1e-10
abs(vy$V_int - VNL_star) < 1e-10

# c) total variance explained equals R2_total
# SHOULD BE TRUE
abs(vx$V_tot - 0.5) < 1e-10
abs(vy$V_tot - 0.5) < 1e-10

# ------------------------------------------------------------------------------------------------------------------
# 4: test that min_abs / max_abs bounds are respected
# ------------------------------------------------------------------------------------------------------------------

# min_abs / max_abs bounds are respected
# SHOULD BE TRUE
all(abs(d$Delta["X", ]) >= 0 & abs(d$Delta["X", ]) <= 1)
all(abs(d$Delta["Y", ]) >= 0 & abs(d$Delta["Y", ]) <= 1)

# stricter example
d_bounds <- sample_delta_1(k = 3,
                           Sigma = Sigma,
                           R2_total = 0.5,
                           R2_interaction = 0.5,
                           include_2way = TRUE,
                           include_3way = TRUE,
                           min_abs = 0.05,
                           max_abs = 1,
                           max_tries = 1000000)

# SHOULD BE TRUE
all(abs(d_bounds$Delta["X", ]) >= 0.05 & abs(d_bounds$Delta["X", ]) <= 1)
all(abs(d_bounds$Delta["Y", ]) >= 0.05 & abs(d_bounds$Delta["Y", ]) <= 1)

# ------------------------------------------------------------------------------------------------------------------
# 5: infeasible inputs should error
# ------------------------------------------------------------------------------------------------------------------

# sampling should run out of tries
sample_delta_1(k = 3,
               Sigma = Sigma,
               R2_total = 0.7,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0.2,
               max_abs = 0.4,
               max_tries = 10000)

# ------------------------------------------------------------------------------------------------------------------
# 6: test that the sampler is reproducible under set.seed()
# ------------------------------------------------------------------------------------------------------------------

# SHOULD BE TRUE
set.seed(1)
d1 <- sample_delta_1(k = 3,
                     Sigma = Sigma,
                     R2_total = 0.5,
                     R2_interaction = 0.5,
                     include_2way = TRUE,
                     include_3way = TRUE,
                     min_abs = 0,
                     max_abs = 1,
                     max_tries = 1000000)

set.seed(1)
d2 <- sample_delta_1(k = 3,
                     Sigma = Sigma,
                     R2_total = 0.5,
                     R2_interaction = 0.5,
                     include_2way = TRUE,
                     include_3way = TRUE,
                     min_abs = 0,
                     max_abs = 1,
                     max_tries = 1000000)

all.equal(d1, d2)
