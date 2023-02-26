# Calculate growth rates (dt) as log ratio to compare with GAM dt

calc_dt_chain <- function(x){
  
  N <- x$N
  
  # calculate population growth rate (chain method)
  dt = c(1) # initial value
  for(i in 2:length(N)){
    dt[i] = log10(N[i]/N[i-1])
  }
  x$dt_chain <- dt
  return(x)
}

# import results of each scenario's LPI
scenarios <- lapply(paste0("simulations/", list.files(path = "simulations/", pattern = "_l.RDS")), readRDS)
names(scenarios) <- gsub("_l.RDS", "", list.files(path = "simulations/", pattern = "_l.RDS"))

# calculate dt via chain
for(i in 1:64){
  
  scenario <- scenarios[[i]] %>% group_by(set, pop) %>% group_split()
  scenario_dt <- lapply(scenario, calc_dt_chain) %>% bind_rows()
  
  saveRDS(scenario_dt, paste0("outputs/", names(scenarios)[i], "_l_dtchain.RDS"))
}
