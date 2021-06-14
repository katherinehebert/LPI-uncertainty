# Script to get the true error introduced into the simulated population time series
# in order to compare index precision to a true expected value

library(dplyr)
simnames <- gsub("_l.RDS", "", list.files(path = "~/Documents/GitHub/LPI-sensitivity/simulations/", pattern = "_l.RDS"))

for(m in 1:length(simnames)){
  
  # get simulation outputs from run_sims.R
  sims <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/simulations/", simnames[m], "_l.RDS")) # simulated pops
  true <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/", simnames[m],"_true.RDS")) # without error
  # add columns to store values later on
  sims$dt <- as.numeric(NA)
  sims$uncertainty <- as.numeric(NA)
  sims$uncertainty_LPI <- as.numeric(NA)
  
  # split into list by population
  sims <- group_by(sims, set, pop) %>% group_split()
  
  # function to calculate the growth rate of the simulated pop time series
  calc_dt <- function(N){
    dt <- 0
    for(n in 2:length(N)){
      # calculate dt
      dt[n] = log10(N[n]/N[n-1])
    }
    return(dt)
  }
  
  # calculate dt for each populations
  for(i in 1:length(sims)){
    sims[[i]]$dt <- calc_dt(sims[[i]]$N)
  }
  
  # combine into one df again to reorganise
  sims <- bind_rows(sims)
  # split into i and j sets because their true trend is differet
  i_sim <- filter(sims, set == "i") %>% group_split(pop)
  j_sim <- filter(sims, set == "j") %>% group_split(pop)
  
  # separate the true sets
  i_true <- filter(true, set == "i")$dt
  j_true <- filter(true, set == "j")$dt
  
  # calculate difference between no-error and with-error growth rates
  # then transform into LPI
  source('~/Documents/GitHub/LPI-sensitivity/scripts/scenario_functions.R')
  for(i in 1:length(i_sim)){ 
    i_sim[[i]]$uncertainty <- i_true - i_sim[[i]]$dt 
    i_sim[[i]]$uncertainty_LPI <- calclpi(i_sim[[i]]$uncertainty)
    }
  for(i in 1:length(j_sim)){
    j_sim[[i]]$uncertainty <- j_true - j_sim[[i]]$dt 
    j_sim[[i]]$uncertainty_LPI <- calclpi(j_sim[[i]]$uncertainty)
    }
  uncertainty <- rbind(bind_rows(i_sim), bind_rows(j_sim))
  
  
  
  # save object
  saveRDS(uncertainty, paste0("outputs/", simnames[m], "_uncertainty.RDS"))
}