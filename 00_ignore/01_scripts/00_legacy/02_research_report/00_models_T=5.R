# In this script we outline the models used to fit the simulated data. 
# note that here we did not include linear control for observed confounders
# as the number of confounders might change. 

# ri-clpm model specification
model_riclpm <- "
rix =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
riy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5
rix ~~ rix
riy ~~ riy
rix ~~ riy

x1 ~~ 0*x1; x2 ~~ 0*x2; x3 ~~ 0*x3; x4 ~~ 0*x4; x5 ~~ 0*x5
y1 ~~ 0*y1; y2 ~~ 0*y2; y3 ~~ 0*y3; y4 ~~ 0*y4; y5 ~~ 0*y5

wx1 =~ 1*x1; wx2 =~ 1*x2; wx3 =~ 1*x3; wx4 =~ 1*x4; wx5 =~ 1*x5
wy1 =~ 1*y1; wy2 =~ 1*y2; wy3 =~ 1*y3; wy4 =~ 1*y4; wy5 =~ 1*y5

rix ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
rix ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5
riy ~~ 0*wx1 + 0*wx2 + 0*wx3 + 0*wx4 + 0*wx5
riy ~~ 0*wy1 + 0*wy2 + 0*wy3 + 0*wy4 + 0*wy5

wx1 ~~ wx1; wx2 ~~ wx2; wx3 ~~ wx3; wx4 ~~ wx4; wx5 ~~ wx5
wy1 ~~ wy1; wy2 ~~ wy2; wy3 ~~ wy3; wy4 ~~ wy4; wy5 ~~ wy5
wy1 ~~ wx1; wy2 ~~ wx2; wy3 ~~ wx3; wy4 ~~ wx4; wy5 ~~ wx5

wx2 ~ wx1 + wy1
wy2 ~ wx1 + wy1
wx3 ~ wx2 + wy2
wy3 ~ wx2 + wy2
wx4 ~ wx3 + wy3
wy4 ~ wx3 + wy3
wx5 ~ wx4 + wy4
wy5 ~ wx4 + wy4

x1 + x2 + x3 + x4 + x5 ~ mx*1
y1 + y2 + y3 + y4 + y5 ~ my*1
"

# dpm model specification
model_dpm <- "
fx =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
fy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5

fx ~~ x1 + y1
fy ~~ x1 + y1

x2 + y2 ~ x1 + y1
x3 + y3 ~ x2 + y2
x4 + y4 ~ x3 + y3
x5 + y5 ~ x4 + y4

x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
x4 ~~ y4
x5 ~~ y5

fx ~~ fx
fy ~~ fy
fx ~~ fy

x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
y1 ~~ y1
y2 ~~ y2
y3 ~~ y3
y4 ~~ y4
y5 ~~ y5

x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
y1 ~ 1
y2 ~ 1
y3 ~ 1
y4 ~ 1
y5 ~ 1
"

# clpm model specification
model_clpm <- "
x2 + y2 ~ x1 + y1
x3 + y3 ~ x2 + y2
x4 + y4 ~ x3 + y3
x5 + y5 ~ x4 + y4

x1 ~~ y1
x2 ~~ y2
x3 ~~ y3
x4 ~~ y4
x5 ~~ y5

x1 ~~ x1
y1 ~~ y1
x2 ~~ x2
y2 ~~ y2
x3 ~~ x3
y3 ~~ y3
x4 ~~ x4
y4 ~~ y4
x5 ~~ x5
y5 ~~ y5

x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
y1 ~ 1
y2 ~ 1
y3 ~ 1
y4 ~ 1
y5 ~ 1
"

# note here that we also want to do BCA sem using
# linear models 
# lasso regression
# random forest regression
# gradient boosting regression
# optional: tabPFN