# test the simulation function

library(dplyr)
library(tidyr)
library(ggplot2)

l2 <- sim_mech(
  n_pairs = 10, timesteps = 4,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = -0.01,
  process = 0, observation = 0,
  K = 100 + 4*c(0:9),
  lag_value = 2,
  "test_lag2"
)

l1 <- sim_mech(
  n_pairs = 10, timesteps = 4,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = -0.01,
  process = 0, observation = 0,
  K = 100 + 4*c(0:9),
  lag_value = 1,
  "test_lag1"
)

l0 <- sim_mech(
  n_pairs = 10, timesteps = 4,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = -0.01,
  process = 0, observation = 0,
  K = 100 + 4*c(0:9),
  lag_value = 0,
  "test_lag0"
)

l1 <- l1 %>% filter(pop == "pop1")
l2 <- l2 %>% filter(pop == "pop1")
l0 <- l0 %>% filter(pop == "pop1")

ggplot(l2) + geom_line(aes(x = time, y = N, col = set)) + 
  theme_linedraw() + labs(title = "2")
ggplot(l1) + geom_line(aes(x = time, y = N, col = set)) + 
  theme_linedraw() + labs(title = "1")
ggplot(l0) + geom_line(aes(x = time, y = N, col = set)) + 
  theme_linedraw() + labs(title = "0")
