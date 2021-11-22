# Function to mechanically simulate two sets of interacting populations and 
# document their properties (plot time series, calculate and plot correlation and covariance)
# to be used as the true population trend

# alpha_ij = 0
# alpha_ji = 0
# K = K_decline
# lag_value = 0
# filename = "test"
# n_pairs = 10 
# timesteps = 10
# N0i = 100
# N0j = 100 
# lambda_i = 1.5
# lambda_j = 1.5

make_true <- function(
  n_pairs = 10, 
  timesteps = 11, 
  N0i = 100, N0j = 100, 
  lambda_i = 1.5, lambda_j = 1.5, 
  alpha_ij, alpha_ji,
  K,
  lag_value,
  simname = filename){
  
  # set seed for randomisations
  set.seed(2)
  
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
  
  # run simulation --------------------------------------------------------------- 
  
  # initialize matrix to store results (population sizes)
  Ni <- as.matrix(rep(N0i, n_pairs))
  Nj <- as.matrix(rep(N0j, n_pairs))
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 1:timesteps){
    
    # population i
    temp_i = Ni[,t]*(1 + lambda_i*(1 - (Ni[,t] + alpha_ij*Nj[,t])/K[t]))

    # if NAs or 0, set to 0.
    temp_i[which(is.na(temp_i))] <- 0
    temp_i[which(temp_i < 0)] <- 0
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
    
    # population j
    temp_j = Nj[,t]*(1 + lambda_j*(1 - (Nj[,t] + alpha_ji*Ni[,t])/K[t]))

    # if NAs or 0, set to 0.
    temp_j[which(is.na(temp_j))] <- 0
    temp_j[which(temp_j < 0)] <- 0
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
  }
  
  # introduce lag to populations j
  if(lag_value != 0){
    # print("running") # to test the loop
    Nj <- Nj[,1:(lag_value + 1)]
    dt_j <- matrix(NA, nrow = n_pairs, ncol = timesteps)
    
    for(t in c(lag_value + 1):timesteps){
      
      # population j
      temp_j = Nj[,t]*(1 + lambda_j*(1 - (Nj[,t] + alpha_ji*Ni[,(t-lag_value)])/K[t]))
      temp_j[which(is.na(temp_j))] <- 0
      temp_j[which(temp_j < 0)] <- 0
      Nj <- cbind(Nj, temp_j)
    }
    
  }

  # remove extra steps from introduced lag
  timesteps <- timesteps - lag_value
  Ni <- Ni[,c(2:(timesteps+1))]
  Nj <- Nj[,c(2:(timesteps+1))]

  # plot results -----------------------------------------------------------------
  
  # create vector of time values for plotting
  time <- 1:timesteps
  
  # function to wrangle the results into long format 
  pops_long <- function(pops_df, n = n_pairs, g = time, set_id) {
    pops_df = as.data.frame(pops_df)
    colnames(pops_df) = as.character(g)
    pops_df$popID <- paste(set_id, sprintf("pop%s", 1:n), sep = "-") 
    pops_df <- pivot_longer(pops_df, 
                            cols = !popID, 
                            names_to = "time", 
                            values_to = "N",
    ) %>% 
      tidyr::separate(popID, into = c("set", "pop"), sep = "-", remove = FALSE) 
    pops_df$time <- as.integer(pops_df$time)
    return(pops_df)
  }
  # bind together
  N = rbind(pops_long(Ni, set_id = "i"), pops_long(Nj, set_id = "j"))
  
  true <- N

  calc_dt_chain <- function(x){
    
    N <- x$N
    
    # calculate population growth rate (chain method)
    dt = c(1) # initial value
    for(i in 2:length(N)){
      dt[i] = log10(N[i]/N[i-1])
    }
    x$dt <- dt
    return(x)
  }
  
  # get annual mean growth rate for each population
  true <- true %>% dplyr::group_by(set, pop) %>% group_split() %>%
    lapply(calc_dt_chain) %>% bind_rows()
  
  # save output
  saveRDS(true, paste0("outputs/", simname, "_true.RDS"))
}