# this script serves to test the delta sampler
#
#
# Recall that the sampler:
# - draws delta coefficients for X and Y at time t = 1,
# - enforces that the sum of squared linear coefficients equals R2_lin,
# - enforces that the sum of squared interaction coefficients equals R2_int,
# - ensures that the total sum of squares equals R2_total,
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
#
# 2) That the returned object has the correct structure:
#    a) A numeric matrix.
#    b) Two rows named "X" and "Y".
#    c) The correct number of columns given k and the toggles.
#    d) Correct column naming and ordering.
#
# 3) That the mathematical constraints are satisfied:
#    a) Sum of squared linear coefficients equals R2_lin.
#    b) Sum of squared interaction coefficients equals R2_int.
#    c) Total sum of squares equals R2_total.
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
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v1_delta_sampler.R"))

# reproducibility
set.seed(1)

# ------------------------------------------------------------------------------------------------------------------
# 1: test that the sampler errors under logically inconsistent inputs
# ------------------------------------------------------------------------------------------------------------------

# a) R2_total outside [0, 1]
# SHOULD ERROR 
sample_delta_1(k = 3,
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
# a) a numeric matrix with 2 rows and k + 1 columns.
# b) two row names "X" and "Y"
# c) the correct number of columns given k and the toggles.
# d) correct column naming and ordering.
d <- sample_delta_1(k = 3,
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

# ------------------------------------------------------------------------------------------------------------------
# 3: test that the mathematical constraints are satisfied
# ------------------------------------------------------------------------------------------------------------------

# a) sum of squared coefficients equals R2_total
# SHOULD BE TRUE
abs(sum(d[1, 1:7]^2) - 0.5) < 1e-10
abs(sum(d[2, 1:7]^2) - 0.5) < 1e-10

# b) sum of squared interaction coefficients equals R2_int
# SHOULD BE TRUE
abs(sum(d[1, 4:7]^2) - 0.25) < 1e-10
abs(sum(d[2, 4:7]^2) - 0.25) < 1e-10

# ------------------------------------------------------------------------------------------------------------------
# 4: test that min_abs / max_abs bounds are respected
# ------------------------------------------------------------------------------------------------------------------

# min_abs / max_abs bounds are respected
# SHOULD BE FALSE
abs(max(d[1, 1:7]) - 1) < 1e-10
abs(min(d[1, 1:7]) - 0) < 1e-10
abs(max(d[2, 1:7]) - 1) < 1e-10
abs(min(d[2, 1:7]) - 0) < 1e-10

# ------------------------------------------------------------------------------------------------------------------
# 5: infeasible inputs should error
# ------------------------------------------------------------------------------------------------------------------

# sampling should run out of tries
sample_delta_1(k = 3,
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
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

set.seed(1)
d2 <- sample_delta_1(k = 3,
               R2_total = 0.5,
               R2_interaction = 0.5,
               include_2way = TRUE,
               include_3way = TRUE,
               min_abs = 0,
               max_abs = 1,
               max_tries = 1000000)

all(d1 == d2)

