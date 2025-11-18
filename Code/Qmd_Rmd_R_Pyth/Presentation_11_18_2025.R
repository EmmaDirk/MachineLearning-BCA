library(tidyverse)
library(patchwork)

# reading data
linear    <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_linear_Pop.csv")
cubic_10  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.10_Pop.csv")
cubic_18  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.18_Pop.csv")
cubic_41  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.41_Pop.csv")
cubic_89  <- read.csv("C:/Users/Admin/Desktop/UU/Thesis/Code/Data/Data_Sim_Cubic_0.89_Pop.csv")

# take samples of each dataset
linear <- linear %>% sample_n(2000)
cubic_10 <- cubic_10 %>% sample_n(2000)
cubic_18 <- cubic_18 %>% sample_n(2000)
cubic_41 <- cubic_41 %>% sample_n(2000)
cubic_89 <- cubic_89 %>% sample_n(2000)

# regress x1 on c1
mod_lin     <- lm(x1 ~ c1, data = linear)
mod_cub_10  <- lm(x1 ~ c1, data = cubic_10)
mod_cub_18  <- lm(x1 ~ c1, data = cubic_18)
mod_cub_41  <- lm(x1 ~ c1, data = cubic_41)
mod_cub_89  <- lm(x1 ~ c1, data = cubic_89)

# add residuals as columns in the data
linear$resid    <- resid(mod_lin)
cubic_10$resid  <- resid(mod_cub_10)
cubic_18$resid  <- resid(mod_cub_18)
cubic_41$resid  <- resid(mod_cub_41)
cubic_89$resid  <- resid(mod_cub_89)

# plot residuals versus predictor
p1 <- ggplot(linear, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Linear")

p2 <- ggplot(cubic_10, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.10")

p3 <- ggplot(cubic_18, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.18")

p4 <- ggplot(cubic_41, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.41")

p5 <- ggplot(cubic_89, aes(x = c1, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal() +
  ggtitle("Cubic 0.89")

# put together with patchwork
(p1 | p2) /
(p3 | p4) /
p5