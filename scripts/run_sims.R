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
max_lambda = 1.5
proc = 0.1
obs = 10

## CARRYING CAPACITY SCENARIOS ##  ------
K_increase = 100 + 4*c(0:9)
K_stable = rep(100, 10)
K_decline = 100 - 4*c(0:9)
# plot these!

## SET A: no covariance from interactions ## ------------------------------------------

## scenario01: extreme decline, no covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0, alpha_ji = 0,
  process = proc, observation = obs,
  K = K_decline,
  lag_value = 0,
  "scenario01"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_decline,
          lag_value = 0,
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
  lag_value = 0,
  "scenario02"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_stable,
          lag_value = 0,
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
  lag_value = 0,
  "scenario03"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0, alpha_ji = 0,
          K = K_increase,
          lag_value = 0,
          "scenario03")
make_gam("scenario03")
get_lpi("scenario03")


## SET B: no-lag covariance from interactions ## -----------------------------------

## scenario 4: extreme decline, negative covariance, no lag  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  lag_value = 0,
  "scenario04"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          lag_value = 0,
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
  lag_value = 0,
  "scenario05"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 0,
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
  lag_value = 0,
  "scenario06"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 0,
          "scenario06")
make_gam("scenario06")
get_lpi("scenario06")


# scenario 7: decline, positive covariance, no lag ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  lag_value = 0,
  "scenario07"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          lag_value = 0,
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
  lag_value = 0,
  "scenario08"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 0,
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
  lag_value = 0,
  "scenario09"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 0,
          "scenario09")
make_gam("scenario09")
get_lpi("scenario09")


## SET C - lag-1 covariance from interactions ----------------------------------

## scenario 10: extreme decline, negative covariance, lag-1  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  lag_value = 1,
  "scenario10"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          lag_value = 0,
          "scenario10")
make_gam("scenario10")
get_lpi("scenario10")


# scenario 11: stable, negative covariance, lag-1 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  lag_value = 1,
  "scenario11"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 1,
          "scenario11")
make_gam("scenario11")
get_lpi("scenario11")

# scenario 12: growth, negative covariance, lag-1 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 1,
  "scenario12"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 1,
          "scenario12")
make_gam("scenario12")
get_lpi("scenario12")


# scenario 13: growth, positive covariance, lag-1 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 1,
  "scenario13"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 1,
          "scenario13")
make_gam("scenario13")
get_lpi("scenario13")

# scenario 14: stable, positive covariance, lag-1 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  lag_value = 1,
  "scenario14"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 1,
          "scenario14")
make_gam("scenario14")
get_lpi("scenario14")

# scenario 15: growth, positive covariance, lag-1 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 1,
  "scenario15"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 1,
          "scenario15")
make_gam("scenario15")
get_lpi("scenario15")


## SET C - lag-2 covariance from interactions ----------------------------------

## scenario 16: extreme decline, negative covariance, lag-2  ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_decline,
  lag_value = 2,
  "scenario16"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_decline,
          lag_value = 0,
          "scenario16")
make_gam("scenario16")
get_lpi("scenario16")


# scenario 17: stable, negative covariance, lag-2 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  lag_value = 2,
  "scenario17"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 2,
          "scenario17")
make_gam("scenario17")
get_lpi("scenario17")

# scenario 18: growth, negative covariance, lag-2 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 2,
  "scenario18"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = -0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 2,
          "scenario18")
make_gam("scenario18")
get_lpi("scenario18")


# scenario 19: growth, positive covariance, lag-2 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 2,
  "scenario19"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 2,
          "scenario19")
make_gam("scenario19")
get_lpi("scenario19")

# scenario 20: stable, positive covariance, lag-2 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_stable,
  lag_value = 2,
  "scenario20"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_stable,
          lag_value = 2,
          "scenario20")
make_gam("scenario20")
get_lpi("scenario20")

# scenario 21: growth, positive covariance, lag-2 ------
sim_mech(
  n_pairs = pop_pairs, timesteps = steps,
  N0i = 100, N0j = 100,
  lambda_i = 1.5, lambda_j = 1.5,
  alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
  process = proc, observation = obs,
  K = K_increase,
  lag_value = 2,
  "scenario21"
)
make_true(n_pairs = pop_pairs, timesteps = steps,
          N0i = 100, N0j = 100,
          lambda_i = 1.5, lambda_j = 1.5,
          alpha_ij = 0.01, alpha_ji = 0.01, # covariance happens here - antagonistic interaction
          K = K_increase,
          lag_value = 2,
          "scenario21")
make_gam("scenario21")
get_lpi("scenario21")