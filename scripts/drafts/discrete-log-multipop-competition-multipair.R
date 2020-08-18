# Multispecies discrete logistic population growth simulation with competition
# for multiple pairs of populations

# simulating growth/decline via changing carrying capacity, but constant growth
# rate and interaction strengths

# load packages
require(tidyverse)
# set ggplot theme
theme_set(theme_linedraw() + theme(panel.grid = element_blank()))

# function to simulate populations

sim_competitors <- function(gen = 10, K0, Kchange){
  
  # gen: number of generations to simulate
  # max_rate: maximum growth rate (r) used for both populations
  
  # constant parameters for simulation -----------------------------------------
  
  # max growth rate
  r_A = 0.1
  r_B = 0.1
  
  # common carrying capacity that changes through time
  K = Kchange*(1:gen) + K0
  K = c(K0, K)
  
  # current population size
  N0_A = 1000
  N0_B = 1000
  
  # interaction coefficients (competition)
  alpha_A = .5
  alpha_B = .9
  
  # run simulation ---------------------------------------------------------------
  
  # initialize vector to store results (population sizes)
  popA <- N0_A
  popB <- N0_B
  
  # initialize variable to hold current population size
  N_A <- N0_A
  N_B <- N0_B
  
  # calculate population sizes
  for(i in 1:gen){
    
    # population A
    N_A <- N_A + N_A*r_A*(1-N_A/K[i] - alpha_A*N_B/K[i])
    if(N_A < 0) {N_A <- 0} # assign 0 to negative population sizes
    popA <- c(popA, N_A) # append resulting population size to results vector
    
    # population B
    N_B <- N_B + N_B*r_B*(1-N_B/K[i] - alpha_B*N_A/K[i])
    if(N_B < 0) {N_B <- 0} # assign 0 to negative population sizes
    popB <- c(popB, N_B) # append resulting population size to results vector
    
  }
  
  # format results as data frame
  results <- data.frame(time = 1:length(popA), K = K, popA = popA, popB = popB)
  
  return(results)
}

# variable simulation parameters

n_pairs = 10 # number of population pairs
K0 = 1000
Kchange = runif(n_pairs, min = -50, max = 0)

sims <- list()
for(i in 1:n_pairs){
  sims[[i]] = sim_competitors(K0 = K0, Kchange = Kchange[i])
}
names(sims) <- sprintf("sim%s",seq(1:n_pairs))

# plot results -----------------------------------------------------------------

# bind results together
sims_df <- bind_rows(sims, .id = "sim") %>%
  pivot_longer(cols = c(popA, popB), names_to = "pop", values_to = "N") %>%
  mutate(id = paste(sim, pop, sep = "_"))
# plot
ggplot(data = sims_df) +
  geom_line(aes(x = time, y = N, group = id, col = pop)) +
  geom_line(aes(x = time, y = K), lty = 2) +
  ylim(min(sims_df$N)-100, max(sims_df$K)+100) +
  facet_wrap(~sim)


# get growth rates -------------------------------------------------------------

getrates <- function(df){
  for(t in 2:nrow(df)){
    df$rA = df$popA[t]/df$popA[t-1] 
    df$rB = df$popB[t]/df$popB[t-1] 
  }
  return(subset(df, select = c(time, K, rA, rB)))
}
# apply to all simulation dataframes
rates = lapply(sims, getrates) %>%
  bind_rows(., .id = "sim") %>%
  pivot_longer(cols = c(rA, rB), names_to = "pop", values_to = "r") %>%
  mutate(id = paste(sim, pop, sep = "_"))
# plot
ggplot(data = rates) +
  geom_line(aes(x = time, y = r, group = id, col = pop)) +
  ylim(min(rates$r)-.05, max(rates$r)+.05) +
  facet_wrap(~sim)
