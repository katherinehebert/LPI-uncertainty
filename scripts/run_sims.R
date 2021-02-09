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

## scenario01: extreme decline, no covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_decline,
  "scenario01"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          "scenario01")
make_gam("scenario01")
get_lpi("scenario01")


# scenario02: stable, no covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_stable,
  "scenario02"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          "scenario02")
make_gam("scenario02")
get_lpi("scenario02")


# scenario 3: growth, no covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_increase,
  "scenario03"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_increase,
          "scenario03")
make_gam("scenario03")
get_lpi("scenario03")


## immediate covariance from interactions ## -----------------------------------

## scenario 4: extreme decline, negative covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  "scenario04"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          "scenario04")
make_gam("scenario04")
get_lpi("scenario04")


# scenario 5: stable, negative covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  "scenario05"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          "scenario05")
make_gam("scenario05")
get_lpi("scenario05")

# scenario 6: growth, negative covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  "scenario06"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          "scenario06")
make_gam("scenario06")
get_lpi("scenario06")


# scenario 7: growth, positive covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  "scenario07"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          "scenario07")
make_gam("scenario07")
get_lpi("scenario07")

# scenario 8: stable, positive covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  "scenario08"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          "scenario08")
make_gam("scenario08")
get_lpi("scenario08")

# scenario 9: growth, positive covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  "scenario09"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          "scenario09")
make_gam("scenario09")
get_lpi("scenario09")