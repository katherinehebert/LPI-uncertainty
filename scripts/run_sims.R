# Script to run simulations for scenarios of population change and calculate LPI 
# for each scenario

# generate the scenarios
source('scripts/scenario_functions.R')
source('scripts/sim_mech.R')
source('scripts/make_true.R')
source('scripts/make_gam.R')
source('scripts/get_lpi.R')

# set theme for all ggplots
theme_set(ggpubr::theme_pubr())

## COMMON PARAMS ## ------
pop_pairs = 10
steps = 10
max_lambda = 2
proc = 0.1
obs = 10

## CARRYING CAPACITY SCENARIOS ##  ------
K_increase = 100 + 4*c(0:9)
K_stable = rep(100, 10)
K_decline = 100 - 4*c(0:9)
# plot these!

## no covariance from interactions ## ------------------------------------------

## scenario 1: extreme decline, no covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_decline,
  "scenario1"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          "scenario1")
make_gam("scenario1")
get_lpi("scenario1")


# scenario 2: stable, no covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_stable,
  "scenario2"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          "scenario2")
make_gam("scenario2")
get_lpi("scenario2")


# scenario 3: growth, no covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_increase,
  "scenario3"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_increase,
          "scenario3")
make_gam("scenario3")
get_lpi("scenario3")


## immediate covariance from interactions ## -----------------------------------

## scenario 4: extreme decline, negative covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  "scenario4"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          "scenario4")
make_gam("scenario4")
get_lpi("scenario4")


# scenario 5: stable, no covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_stable,
  "scenario2"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          "scenario2")
make_gam("scenario2")
get_lpi("scenario2")

# 
# # scenario 3: growth, no covariance, no lag ------
# sim_mech(
#   n_pairs = pop_pairs, timesteps = steps,
#   N0i = 100, N0j = 100,
#   lambda_i = 1.5, lambda_j = 1.5,
#   alpha_ij = 0, alpha_ji = 0,
#   process = proc, observation = obs,
#   K = K_increase,
#   "scenario3"
# )
# make_true(n_pairs = pop_pairs, timesteps = steps,
#           N0i = 100, N0j = 100,
#           lambda_i = 1.5, lambda_j = 1.5,
#           alpha_ij = 0, alpha_ji = 0,
#           K = K_increase,
#           "scenario3")
# make_gam("scenario3")
# get_lpi("scenario3")
