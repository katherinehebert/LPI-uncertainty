# Function to smechanically simulate two sets of interacting populations and 
# document their properties (plot time series, calculate and plot correlation and covariance)

mech_sim <- function(
  npairs, 
  timesteps, 
  N0i = 200, N0j = 200, 
  lambda_i, lambda_j, 
  alpha_ij, alpha_ji,
  process, 
  observation,
  K,
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
  # process = absolute value of the process error limit to add to growth rates
  # observation = absolute value of the obs error limit to add to population growth rates
  
  # make growth rates matrix form
  r_i = matrix(lambda_i, nrow = n_pairs, ncol = timesteps) 
  r_j = matrix(lambda_i, nrow = n_pairs, ncol = timesteps) 
  
  # add process error to growth rates
  process = rnorm(n = n_pairs*timesteps, mean = process, sd = process/10)
  r_i_error = r_i + process
  r_j_error = r_j + process
  # observation error
  obs_i_error = matrix(rnorm(n = n_pairs*timesteps, mean = observation, sd = observation/10), nrow = n_pairs, ncol = timesteps)
  obs_j_error = matrix(rnorm(n = n_pairs*timesteps, mean = observation, sd = observation/10), nrow = n_pairs, ncol = timesteps)
  
  # run simulation --------------------------------------------------------------- 
  
  # initialize matrix to store results (population sizes)
  Ni <- as.matrix(rep(N0i, n_pairs))
  Nj <- as.matrix(rep(N0j, n_pairs))
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 1:timesteps-1){
    
    # population i
    temp_i = Ni[t]*(1 + r_i_error[,t]*(1 - (Ni[t] + alpha_ij*Nj[t])/K)) + obs_i_error[,t]
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
    
    # population j
    temp_j = Nj[t]*(1 + r_j_error[,t]*(1 - (Nj[t] + alpha_ji*Ni[t])/K)) + obs_j_error[,t]
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
  }
  
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
  
  # plot
  ggplot(N) +
    geom_line(aes(x = time, y = N, group = popID, col = popID)) + 
    #geom_hline(yintercept = K, lty = 2) +
    #geom_hline(yintercept = K, lty = 2) +
    facet_wrap(~ set) + 
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max(N$N)+10))
  
  # save outputs -----------------------------------------------------------------
  saveRDS(N, paste0("simulations/", simname, "_l.RDS"))
  ggsave(filename = paste0(simname, "_N.png"), path = "figures/", plot = last_plot(),
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
  
  
  # # make parameter table
  # params <- data.frame("i" = c(N0i, 
  #                              r_i[1,1], 
  #                              alpha_ji, 
  #                              paste0("+/-", observation), 
  #                              paste0("+/-", process)), 
  #                      "j" = c(N0j, 
  #                              r_j[1,1], 
  #                              alpha_ji, 
  #                              paste0("+/-", observation), 
  #                              paste0("+/-", process)),
  #                      row.names = c("Initial size (N0)",
  #                                    "Maximum growth rate (r)",
  #                                    "Interaction effect (\\alpha)",
  #                                    "Observation error",
  #                                    "Process error"))
  # write.table(params, paste0("simulations/", simname, "_params.txt"))
  # 
  return(N)
}
