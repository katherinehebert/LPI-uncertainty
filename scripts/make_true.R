# Function to mechanically simulate two sets of interacting populations and 
# document their properties (plot time series, calculate and plot correlation and covariance)
# to be used as the true population trend

make_true <- function(
  n_pairs, 
  timesteps, 
  N0i, N0j, 
  lambda_i, lambda_j, 
  alpha_ij, alpha_ji,
  K,
  lag_value,
  simname){
  
  ## arguments: ## ----
  
  # npairs = number of pairs of populations to simulate
  # timesteps = number of time steps
  # N0i = initial population size for set i
  # N0j = initial population size for set j
  # lambda_i = true population growth rate for set i
  # lambda_j = true population growth rate for set j
  # alpha_ij = interaction coefficient (effect of set j on set i)
  # alpha_ji = interaction coefficient (effect of set i on set j)
  # K = carrying capacity of the environment
  
  # add more time steps to allow lag
  timesteps = timesteps + lag_value
  
  # make growth rates matrix form
  r_i = matrix(lambda_i, nrow = n_pairs, ncol = timesteps) 
  r_j = matrix(lambda_i, nrow = n_pairs, ncol = timesteps) 
  
  # run simulation --------------------------------------------------------------- 
  
  # initialize matrix to store results (population sizes)
  Ni <- as.matrix(rep(N0i, n_pairs))
  Nj <- as.matrix(rep(N0j, n_pairs))
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 1:timesteps-1){
    
    t_lag = t - lag_value
    
    # population i
    temp_i = Ni[t]*(1 + r_i[,t]*(1 - (Ni[t] + alpha_ij*Nj[t_lag])/K[t])) 
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
    
    # population j
    temp_j = Nj[t]*(1 + r_j[,t]*(1 - (Nj[t] + alpha_ji*Ni[t])/K[t])) 
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
  }
  
  # remove extra steps from introduced lag
  timesteps = timesteps-lag_value
  Ni <- Ni[,c(1:timesteps)]
  Nj <- Nj[,c(1:timesteps)]
  
  # plot results -----------------------------------------------------------------
  
  # create vector of time values for plotting
  time <- 1:timesteps
  
  # function to wrangle the results into long format 
  pops_long <- function(pops_df, n = n_pairs, g = timesteps, set_id) {
    pops_df = as.data.frame(pops_df)
    colnames(pops_df) = time
    pops_df = mutate(.data = pops_df, "popID" = paste(set_id, sprintf("pop%s", 1:n), sep = "-")) %>%
      pivot_longer(cols = 1:all_of(g), names_to = "time", values_to = "N") %>%
      separate(popID, into = c("set", "pop"), sep = "-", remove = FALSE) %>%
      mutate_at(vars(time), as.integer)
  }
  # bind together
  N = rbind(pops_long(Ni, set_id = "i"), pops_long(Nj, set_id = "j"))
  
  true <- filter(N, pop == "pop1")
  true$dt <- NA  
  for(i in 2:10){
    x <- filter(true, popID == "i-pop1")
    y <- filter(true, popID == "j-pop1")
    true$dt[i] <- log10(x$N[i]/x$N[i-1])
    true$dt[i+10] <- log10(y$N[i]/y$N[i-1])
  }
  saveRDS(true, paste0("outputs/", simname, "_true.RDS"))
}