# Function to mechanically simulate two sets of interacting populations and 
# document their properties (plot time series, calculate and plot correlation and covariance)

library(dplyr)
library(tidyr)
library(ggplot2)

sim_mech <- function(
  n_pairs = 10, 
  timesteps = 10, 
  lambda_i = max_lambda, lambda_j = max_lambda, 
  alpha_ij = 0, alpha_ji = 0,
  process = proc, 
  observation = obs,
  K,
  lag_value,
  simname = filename,
  save_figs = TRUE){
  
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
  # process = absolute value of the  to add to growth rates
  # observation = absolute value of the obs error limit to add to population sizes
  
  # add more time steps to allow lag
  timesteps = timesteps + lag_value
  
  # calculate initial population size as the equilibrium population size -------
  
  # calculate starting population size at equilibrium
  get_equilibrium = function(K0, alpha_1, alpha_2){
    N_eq = (K0 - alpha_1 * K0)/(1 - alpha_2 * alpha_1)
    return(N_eq)
  }
  N0i = get_equilibrium(K0 = K[1],
                        alpha_1 = alpha_ij,
                        alpha_2 = alpha_ji) %>% ceiling()
  N0j = get_equilibrium(K0 = K[1],
                        alpha_1 = alpha_ji,
                        alpha_2 = alpha_ij) %>% ceiling()
  
  # run simulation --------------------------------------------------------------- 
  
  # initialize matrix to store results (population sizes)
  Ni <- as.matrix(rep(N0i, n_pairs))
  Nj <- as.matrix(rep(N0j, n_pairs))
  
  # initialize matrix to store growth rates
  dt_i <- matrix(NA, nrow = n_pairs, ncol = timesteps)
  dt_j <- matrix(NA, nrow = n_pairs, ncol = timesteps)
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 1:timesteps){
    
    # generate growth rates from a lognormal distribution with process error
    r_i_error = rlnorm(n = n_pairs, meanlog = log(lambda_i), sdlog = process)
    r_j_error = rlnorm(n = n_pairs, meanlog = log(lambda_j), sdlog = process)
    
    # population i
    temp_i = Ni[,t]*(1 + r_i_error*(1 - (Ni[,t] + alpha_ij*Nj[,t])/K[t]))

    # if NAs or 0, set to 0.
    temp_i[which(is.na(temp_i))] <- 0
    temp_i[which(temp_i < 0)] <- 0
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
     
    # population j
    temp_j = Nj[,t]*(1 + r_j_error*(1 - (Nj[,t] + alpha_ji*Ni[,t])/K[t]))

    # if NAs or 0, set to 0.
    temp_j[which(is.na(temp_j))] <- 0
    temp_j[which(temp_j < 0)] <- 0
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector

    ## save generated growth rates
    dt_i[,t] <- r_i_error
    dt_j[,t] <- r_j_error
  }
  
  # introduce lag to populations j
  if(lag_value != 0){
    # print("running") # to test the loop
    Nj <- Nj[,1:(lag_value + 1)]

    for(t in c(lag_value + 1):timesteps){
      
      # population j
      temp_j = Nj[,t]*(1 + dt_j[t]*(1 - (Nj[,t] + alpha_ji*Ni[,(t-lag_value)])/K[t]))
      temp_j[which(is.na(temp_j))] <- 0
      temp_j[which(temp_j < 0)] <- 0
      Nj <- cbind(Nj, temp_j)
    }

  }
  
  # apply observation error on the calculated population sizes
  # pick number of lognormal distribution with a mean of N and an sd of observation error
  # but ensure no zeros are introduced from measurement error
  nozeros = function(x) {
    if(x <= 0){x = 1} # correct because cannot log 0 below
    xerr = rlnorm(1, meanlog = log(x), sdlog = observation)
    if(xerr <= 0){
      # keep randomizing until x is no longer zero or below
      while(xerr <= 0){
        xerr = rlnorm(1, meanlog = log(x), sdlog = observation)
      }
    }
    return(xerr)
  }
  Ni <- apply(Ni, 1:2, nozeros) 
  Nj <- apply(Nj, 1:2, nozeros) 

  # remove extra steps from introduced lag
  timesteps <- timesteps - lag_value
  Ni <- Ni[,c(2:(timesteps+1))]
  Nj <- Nj[,c(2:(timesteps+1))]
  dt_j <- dt_j[,c(1:timesteps)]
  
  # save growth rates
  dt <- list("dt_i" = dt_i, "dt_j" = dt_j)
  saveRDS(dt, paste0("simulations/", simname, "_rerror.RDS"))
  
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
  
  # make parameter table
  params <- data.frame("parameter" = c("Initial size (N0)",
                                       "Maximum growth rate (r)",
                                       "Interaction effect (alpha)",
                                       "Observation error",
                                       "Process error",
                                       "Lag"),
                       "i" = c(N0i,
                               lambda_i,
                               alpha_ij,
                               observation,
                               process,
                               lag_value),
                       "j" = c(N0j,
                               lambda_j,
                               alpha_ji,
                               observation,
                               process,
                               lag_value)
  )
  saveRDS(params, paste0("simulations/", simname, "_params.RDS"))
  saveRDS(N, paste0("simulations/", simname, "_l.RDS"))
  
  if(save_figs == FALSE){ 
    return(N)
  } else {
    
    # plot
    N_plot <- ggplot(N) +
      geom_line(aes(x = time, y = N, group = popID, col = popID)) + 
      facet_wrap(~ set) + 
      theme(legend.position = "none") +
      scale_x_continuous(breaks = c(1:10)) + 
      coord_cartesian(ylim = c(0, max(N$N)+10))
    
    # save outputs -----------------------------------------------------------------
    ggsave(filename = paste0(simname, "_N.png"), path = "figures/", plot = N_plot,
           width = 7, height = 5, units = "in")

  
  # calculate covariation --------------------------------------------------------
  
  N_w = cbind(t(Ni), t(Nj))
  
  # plot covariation
  png(paste0("figures/", simname, "_cov.png"), width = 500, height = 500)
  cov(N_w) %>% heatmap(Colv = NA, Rowv = NA, 
                       col = (colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)),
                       main = "Covariation")
  dev.off()
  
  # plot correlation
  png(paste0("figures/", simname, "_cor.png"), width = 500, height = 500)
  cor(N_w) %>% heatmap(Colv = NA, Rowv = NA, 
                       col = (colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)), 
                       main = "Correlation")
  dev.off()
  
  return(N)
    }
}
