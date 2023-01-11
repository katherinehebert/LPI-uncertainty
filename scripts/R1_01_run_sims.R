# Script to run simulations for R1_scenarios of population change and calculate LPI 
# for each R1_scenario

# generate the R1_scenarios
source('scripts/R1_scenario_functions.R')
source('scripts/sim_mech.R')
source('scripts/make_true.R')
source('scripts/make_gam.R')
source('scripts/get_lpi.R')

# set theme for all ggplots
theme_set(ggpubr::theme_pubr())

## COMMON PARAMS ---------------------------------------------------------------
pop_pairs = 1000
steps = 11
obs = 5

## CARRYING CAPACITY R1_scenarioS -------------------------------------------------
K_increase = 100 + 10*c(0:10)
K_stable = rep(100, 11)
K_decline = 100 - 5*c(0:10)
# plot these!
K_plot <- data.frame(
  time = 1:steps,
  R1_scenario = factor(c(rep("increase", steps), 
                      rep("stable", steps), 
                      rep("decline", steps)),
                    levels = rev(c("increase", "stable", "decline"))),
  K = c(K_increase, K_stable, K_decline)
)
ggplot(K_plot) +
  geom_line(aes(x = time, y = K, col = R1_scenario)) +
  labs(x = "", y = "Carrying capacity (K)", 
       col = "Trend") +
  scale_color_manual(values = ggsci::pal_locuszoom("default")(6)[c(1,5,3)]) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1:10))
ggsave("figures/carryingcapacity.png", width = 2.63, height = 2.38)


## SET 1: no covariance from interactions --------------------------------------

# how does direction affect the LPI?
# R1_scenario1A: extreme decline, no covariance, no lag, no process error
# R1_scenario1B: stable, no covariance, no lag, no process error
# R1_scenario1C: growth, no covariance, no lag, no process error 

# is this finding robust to process error?
# level 1: low process error
# R1_scenario1D: extreme decline, no covariance, no lag, low process error
# R1_scenario1E: stable, no covariance, no lag, low process error
# R1_scenario1F: growth, no covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario1G: extreme decline, no covariance, no lag, high process error
# R1_scenario1H: stable, no covariance, no lag, high process error
# R1_scenario1I: growth, no covariance, no lag, high process error

sim_names <- paste0("R1_scenario1", LETTERS[1:9])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 3)
proc_error <- rep(c(0, 0.1, 0.2), each = 3)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = 0, alpha_ji = 0,
    process = proc_error[i],
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
}

## SET 2: negative covariance from interactions --------------------------------

# how does negative covariation affect the LPI?
# R1_scenario2A: decline, negative covariance, no lag, no process error
# R1_scenario2B: stable, negative covariance, no lag, no process error
# R1_scenario2C: growth, negative covariance, no lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario2D: decline, negative covariance, no lag, low process error
# R1_scenario2E: stable, negative covariance, no lag, low process error
# R1_scenario2F: growth, negative covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario2G: decline, negative covariance, no lag, high process error
# R1_scenario2H: stable, negative covariance, no lag, high process error
# R1_scenario2I: growth, negative covariance, no lag, high process error

# R1_scenarios 2J-2R: same levels, but with strong negative covariance

sim_names <- paste0("R1_scenario2", LETTERS[1:18])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
alphas <- rep(list(c(-0.1, 0.1), c(-0.2, 0.2)), each = 9)
proc_error <- rep(rep(c(0, 0.1, 0.2), each = 3), 2)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i],
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
}


## SET 3: positive covariance from interactions --------------------------------

# how does positive covariation affect the LPI?
# R1_scenario3A: extreme decline, positive covariance, no lag, no process error
# R1_scenario3B: stable, positive covariance, no lag, no process error
# R1_scenario3C: growth, positive covariance, no lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario3D: extreme decline, positive covariance, no lag, low process error
# R1_scenario3E: stable, positive covariance, no lag, low process error
# R1_scenario3F: growth, positive covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario3G: extreme decline, positive covariance, no lag, high process error
# R1_scenario3H: stable, positive covariance, no lag, high process error
# R1_scenario3I: growth, positive covariance, no lag, high process error

# R1_scenarios 3J-3R: same levels, but with strong positive covariance

sim_names <- paste0("R1_scenario3", LETTERS[1:18])
alphas <- rep(list(c(-0.1, -0.1), c(-0.2, -0.2)), each = 9)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1],
    alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i],
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
}

## SET 4: LAG-1 negative covariance from interactions --------------------------

# how does lagged negative covariation affect the LPI?
# R1_scenario4A: extreme decline, negative covariance, 1 step lag, no process error
# R1_scenario4B: stable, negative covariance, 1 step lag, no process error
# R1_scenario4C: growth, negative covariance, 1 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario4D: extreme decline, negative covariance, 1 step lag, low process error
# R1_scenario4E: stable, negative covariance, 1 step lag, low process error
# R1_scenario4F: growth, negative covariance, 1 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario4G: extreme decline, negative covariance, 1 step lag, high process error
# R1_scenario4H: stable, negative covariance, 1 step lag, high process error
# R1_scenario4I: growth, negative covariance, 1 step lag, high process error

# R1_scenarios 4J-R: same levels, but with strong negative covariance

sim_names <- paste0("R1_scenario4", LETTERS[1:18])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
alphas <- rep(list(c(-0.1, 0.1), c(-0.2, 0.2)), each = 9)
proc_error <- rep(rep(c(0, 0.1, 0.2), each = 3), 2)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i], 
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 1
  )
}


## SET 5: LAG-1 positive covariance from interactions --------------------------

# how does lagged positive covariation affect the LPI?
# R1_scenario5A: extreme decline, positive covariance, 1 step lag, no process error
# R1_scenario5B: stable, positive covariance, 1 step lag, no process error
# R1_scenario5C: growth, positive covariance, 1 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario5D: extreme decline, positive covariance, 1 step lag, low process error
# R1_scenario5E: stable, positive covariance, 1 step lag, low process error
# R1_scenario5F: growth, positive covariance, 1 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario5G: extreme decline, positive covariance, 1 step lag, high process error
# R1_scenario5H: stable, positive covariance, 1 step lag, high process error
# R1_scenario5I: growth, positive covariance, 1 step lag, high process error

# R1_scenarios 5J-R: same levels, but with strong positive covariance

sim_names <- paste0("R1_scenario5", LETTERS[1:18])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
alphas <- rep(list(c(-0.1, -0.1), c(-0.2, -0.2)), each = 9)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i], 
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 1
  )
}

## SET 6: LAG-2 negative covariance from interactions --------------------------

# how does lagged negative covariation affect the LPI?
# R1_scenario6A: extreme decline, negative covariance, 2 step lag, no process error
# R1_scenario6B: stable, negative covariance, 2 step lag, no process error
# R1_scenario6C: growth, negative covariance, 2 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario6D: extreme decline, negative covariance, 2 step lag, low process error
# R1_scenario6E: stable, negative covariance, 2 step lag, low process error
# R1_scenario6F: growth, negative covariance, 2 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario6G: extreme decline, negative covariance, 2 step lag, high process error
# R1_scenario6H: stable, negative covariance, 2 step lag, high process error
# R1_scenario6I: growth, negative covariance, 2 step lag, high process error

# R1_scenarios 6J-R: same levels, but with strong negative covariance

sim_names <- paste0("R1_scenario6", LETTERS[1:18])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
alphas <- rep(list(c(-0.1, 0.1), c(-0.2, 0.2)), each = 9)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i], 
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 2
  )
}


## SET 7: LAG-2 positive covariance from interactions ## -----------------------------------

# how does lagged positive covariation affect the LPI?
# R1_scenario7A: extreme decline, positive covariance, 2 step lag, no process error
# R1_scenario7B: stable, positive covariance, 2 step lag, no process error
# R1_scenario7C: growth, positive covariance, 2 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# R1_scenario7D: extreme decline, positive covariance, 2 step lag, low process error
# R1_scenario7E: stable, positive covariance, 2 step lag, low process error
# R1_scenario7F: growth, positive covariance, 2 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# R1_scenario7G: extreme decline, positive covariance, 2 step lag, high process error
# R1_scenario7H: stable, positive covariance, 2 step lag, high process error
# R1_scenario7I: growth, positive covariance, 2 step lag, high process error

# R1_scenarios 7J-R: same levels, but with strong positive covariance

sim_names <- paste0("R1_scenario7", LETTERS[1:18])
K_R1_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
alphas <- rep(list(c(-0.1, -0.1), c(-0.2, -0.2)), each = 9)

for(i in 1:length(sim_names)){
  filename <- sim_names[[i]]
  sim_mech(
    n_pairs = pop_pairs, timesteps = steps,
    N0i = 100, N0j = 100,
    lambda_i = 1.5, lambda_j = 1.5,
    alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
    process = proc_error[i], 
    observation = obs,
    K = K_R1_scenarios[[i]],
    lag_value = 2
  )
}
