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

## COMMON PARAMS ---------------------------------------------------------------
pop_pairs = 1000
steps = 11
obs = 5

## CARRYING CAPACITY SCENARIOS -------------------------------------------------
K_increase = 100 + 10*c(0:10)
K_stable = rep(100, 11)
K_decline = 100 - 5*c(0:10)
# plot these!
K_plot <- data.frame(
  time = 1:steps,
  scenario = factor(c(rep("increase", steps), 
                      rep("stable", steps), 
                      rep("decline", steps)),
                    levels = rev(c("increase", "stable", "decline"))),
  K = c(K_increase, K_stable, K_decline)
)
ggplot(K_plot) +
  geom_line(aes(x = time, y = K, col = scenario)) +
  labs(x = "", y = "Carrying capacity (K)", 
       col = "Trend") +
  scale_color_manual(values = ggsci::pal_locuszoom("default")(6)[c(1,5,3)]) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(1:10))
ggsave("figures/carryingcapacity.png", width = 2.63, height = 2.38)


## SET 1: no covariance from interactions --------------------------------------

# how does direction affect the LPI?
# scenario1A: extreme decline, no covariance, no lag, no process error
# scenario1B: stable, no covariance, no lag, no process error
# scenario1C: growth, no covariance, no lag, no process error 

# is this finding robust to process error?
# level 1: low process error
# scenario1D: extreme decline, no covariance, no lag, low process error
# scenario1E: stable, no covariance, no lag, low process error
# scenario1F: growth, no covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario1G: extreme decline, no covariance, no lag, high process error
# scenario1H: stable, no covariance, no lag, high process error
# scenario1I: growth, no covariance, no lag, high process error

sim_names <- paste0("scenario1", LETTERS[1:9])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 3)
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
    K = K_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = 0, alpha_ji = 0,
            process = proc_error[i],
            K = K_scenarios[[i]],
            lag_value = 0)
  # make_gam(filename)
  # get_lpi(filename)
}

## SET 2: negative covariance from interactions --------------------------------

# how does negative covariation affect the LPI?
# scenario2A: decline, negative covariance, no lag, no process error
# scenario2B: stable, negative covariance, no lag, no process error
# scenario2C: growth, negative covariance, no lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario2D: decline, negative covariance, no lag, low process error
# scenario2E: stable, negative covariance, no lag, low process error
# scenario2F: growth, negative covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario2G: decline, negative covariance, no lag, high process error
# scenario2H: stable, negative covariance, no lag, high process error
# scenario2I: growth, negative covariance, no lag, high process error

# scenarios 2J-2R: same levels, but with strong negative covariance

sim_names <- paste0("scenario2", LETTERS[1:18])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
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
    K = K_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            process = proc_error[i],
            K = K_scenarios[[i]],
            lag_value = 0)
  # make_gam(filename)
  # get_lpi(filename)
}


## SET 3: positive covariance from interactions --------------------------------

# how does positive covariation affect the LPI?
# scenario3A: extreme decline, positive covariance, no lag, no process error
# scenario3B: stable, positive covariance, no lag, no process error
# scenario3C: growth, positive covariance, no lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario3D: extreme decline, positive covariance, no lag, low process error
# scenario3E: stable, positive covariance, no lag, low process error
# scenario3F: growth, positive covariance, no lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario3G: extreme decline, positive covariance, no lag, high process error
# scenario3H: stable, positive covariance, no lag, high process error
# scenario3I: growth, positive covariance, no lag, high process error

# scenarios 3J-3R: same levels, but with strong positive covariance

sim_names <- paste0("scenario3", LETTERS[1:18])
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
    K = K_scenarios[[i]],
    lag_value = 0,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1],
            alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            K = K_scenarios[[i]],
            process = proc_error[i],
            lag_value = 0)
  # make_gam(filename)
  # get_lpi(filename)
}

## SET 4: LAG-1 negative covariance from interactions --------------------------

# how does lagged negative covariation affect the LPI?
# scenario4A: extreme decline, negative covariance, 1 step lag, no process error
# scenario4B: stable, negative covariance, 1 step lag, no process error
# scenario4C: growth, negative covariance, 1 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario4D: extreme decline, negative covariance, 1 step lag, low process error
# scenario4E: stable, negative covariance, 1 step lag, low process error
# scenario4F: growth, negative covariance, 1 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario4G: extreme decline, negative covariance, 1 step lag, high process error
# scenario4H: stable, negative covariance, 1 step lag, high process error
# scenario4I: growth, negative covariance, 1 step lag, high process error

# scenarios 4J-R: same levels, but with strong negative covariance

sim_names <- paste0("scenario4", LETTERS[1:18])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
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
    K = K_scenarios[[i]],
    lag_value = 1,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            K = K_scenarios[[i]],
            process = proc_error[i],
            lag_value = 1)
  # make_gam(filename)
  # get_lpi(filename)
}


## SET 5: LAG-1 positive covariance from interactions --------------------------

# how does lagged positive covariation affect the LPI?
# scenario5A: extreme decline, positive covariance, 1 step lag, no process error
# scenario5B: stable, positive covariance, 1 step lag, no process error
# scenario5C: growth, positive covariance, 1 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario5D: extreme decline, positive covariance, 1 step lag, low process error
# scenario5E: stable, positive covariance, 1 step lag, low process error
# scenario5F: growth, positive covariance, 1 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario5G: extreme decline, positive covariance, 1 step lag, high process error
# scenario5H: stable, positive covariance, 1 step lag, high process error
# scenario5I: growth, positive covariance, 1 step lag, high process error

# scenarios 5J-R: same levels, but with strong positive covariance

sim_names <- paste0("scenario5", LETTERS[1:18])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
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
    K = K_scenarios[[i]],
    lag_value = 1,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            process = proc_error[i],
            K = K_scenarios[[i]],
            lag_value = 1)
  # make_gam(filename)
  # get_lpi(filename)
}

## SET 6: LAG-2 negative covariance from interactions --------------------------

# how does lagged negative covariation affect the LPI?
# scenario6A: extreme decline, negative covariance, 2 step lag, no process error
# scenario6B: stable, negative covariance, 2 step lag, no process error
# scenario6C: growth, negative covariance, 2 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario6D: extreme decline, negative covariance, 2 step lag, low process error
# scenario6E: stable, negative covariance, 2 step lag, low process error
# scenario6F: growth, negative covariance, 2 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario6G: extreme decline, negative covariance, 2 step lag, high process error
# scenario6H: stable, negative covariance, 2 step lag, high process error
# scenario6I: growth, negative covariance, 2 step lag, high process error

# scenarios 6J-R: same levels, but with strong negative covariance

sim_names <- paste0("scenario6", LETTERS[1:18])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
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
    K = K_scenarios[[i]],
    lag_value = 2,
    save_figs = FALSE
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            K = K_scenarios[[i]],
            process = proc_error[i],
            lag_value = 2)
  # make_gam(filename)
  # get_lpi(filename)
}


## SET 7: LAG-2 positive covariance from interactions ## -----------------------------------

# how does lagged positive covariation affect the LPI?
# scenario7A: extreme decline, positive covariance, 2 step lag, no process error
# scenario7B: stable, positive covariance, 2 step lag, no process error
# scenario7C: growth, positive covariance, 2 step lag, no process error

# is this finding robust to process error?
# level 1: low process error
# scenario7D: extreme decline, positive covariance, 2 step lag, low process error
# scenario7E: stable, positive covariance, 2 step lag, low process error
# scenario7F: growth, positive covariance, 2 step lag, low process error

# is this finding robust to process error?
# level 2: high process error
# scenario7G: extreme decline, positive covariance, 2 step lag, high process error
# scenario7H: stable, positive covariance, 2 step lag, high process error
# scenario7I: growth, positive covariance, 2 step lag, high process error

# scenarios 7J-R: same levels, but with strong positive covariance

sim_names <- paste0("scenario7", LETTERS[1:18])
K_scenarios <- rep(list(K_decline, K_stable, K_increase), 6)
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
    K = K_scenarios[[i]],
    lag_value = 2,
    save_figs = FALSE 
  )
  make_true(n_pairs = pop_pairs, timesteps = steps,
            N0i = 100, N0j = 100,
            lambda_i = 1.5, lambda_j = 1.5,
            alpha_ij = alphas[[i]][1], alpha_ji = alphas[[i]][2], # this is where the covariance is introduced
            K = K_scenarios[[i]],
            process = proc_error[i],
            lag_value = 2)
  # make_gam(filename)
  # get_lpi(filename)
}
