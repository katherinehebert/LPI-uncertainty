# Functions needed to run the LPI sensitivity tests in the 
# scenario1_predation and scenario2_competition documents.


## GENERAL USE FUNCTIONS ## ----------------------------------------------------

# function to repeat vector as rows of a matrix
rep.row <- function(x,n){matrix(rep(x, each = n), nrow = n)}


## SIMULATION FUNCTIONS ## -----------------------------------------------------

# function to wrangle simulation results into long format 
pops_long <- function(pops_df, n = n_pairs, g = steps, set_id) {
  pops_df = as.data.frame(pops_df)
  colnames(pops_df) = as.character(1:steps)
  pops_df = mutate(.data = pops_df, "popID" = paste(set_id, sprintf("pop%s", 1:n), sep = "-")) %>%
    pivot_longer(cols = 1:all_of(g), names_to = "time", values_to = "N") %>%
    separate(popID, into = c("set", "pop"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(time), as.integer)
}

# function to simulate populations according to Lotka-Volterra competition model
sim <- function(Ni_init = N0i, Nj_init = N0j, n_pairedpops = n_pairs, timesteps = steps, Ri = r_i_error, Rj = r_j_error, int_ij = alpha_ij, int_ji = alpha_ji, proc = proc_error, K){
  
  # initialize matrix to store results (population sizes)
  Ni <- as.matrix(rep(Ni_init, n_pairedpops))
  Nj <- as.matrix(rep(Nj_init, n_pairedpops))
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 2:timesteps-1){
    
    # population i
    temp_i = Ni[t]*(1 + Ri[,t]*(1 - (Ni[t] + int_ij*Nj[t])/K[t])) + proc[,t]
    temp_i[which(temp_i < 0)] <- 0 # assign 0 to negative population sizes
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
    
    # population j
    temp_j = Nj[t]*(1 + Rj[,t]*(1 - (Nj[t] + int_ji*Ni[t])/K[t])) + proc[,t]
    temp_j[which(temp_j < 0)] <- 0 # assign 0 to negative population sizes
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
  }
  
  # wrangle!
  time <- 1:timesteps   # create vector of time values for plotting
  # bind together
  N = rbind(pops_long(Ni, set_id = "i"), pops_long(Nj, set_id = "j"))
  
  return(N)
}

# plot simulated populations
sim_plot <- function(df, K_slope, K_int, title){
  ggplot(df) +
    geom_line(aes(x = time, y = N, group = popID, col = popID)) + 
    geom_abline(slope = K_slope, intercept = K_int, lty = 2, lwd = .6) + 
    facet_wrap(~ set) + 
    ylim(c(0, K_int*2)) +
    ggtitle(title)
}

# calculate and plot covariation and correlation heatmaps
covcor_pops <- function(df_long){
  # convert to wide
  temp_w <- subset(df_long, select = -c(set, pop)) %>%
    pivot_wider(names_from = popID, values_from = N, id_cols = time) %>%
    subset(select = -time)
  # calculate and plot covariance
  cov_res <- cov(temp_w) %>% as.data.frame() %>% mutate(population_j = rownames(.)) %>%
    pivot_longer(cols = 1:(ncol(.)-1), names_to = "population_i") 
  # keep rows of i vs. j only (not i vs. i or j vs. j)
  cov_res = cov_res[-grep("i-", cov_res$population_j),]
  cov_res = cov_res[-grep("j-", cov_res$population_i),]
  # plot heatmap
  p1 <- ggplot(data = cov_res, aes(x = population_j, y = population_i)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis_c() +
    labs(y = "", x = "", fill = "Covariation") + theme(legend.position = "right")
  
  # calculate and plot correlation
  cor_res <- cor(temp_w) %>% as.data.frame() %>% mutate(population_j = rownames(.)) %>%
    pivot_longer(cols = 1:(ncol(.)-1), names_to = "population_i")
  # keep rows of i vs. j only (not i vs. i or j vs. j)
  cor_res = cor_res[-grep("i-", cor_res$population_j),]
  cor_res = cor_res[-grep("j-", cor_res$population_i),]
  
  # plot heatmap
  p2 <- ggplot(data = cor_res, aes(x = population_j, y = population_i)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis_c() + 
    labs(y = "", x = "", fill = "Correlation") + theme(legend.position = "right")
  return(list(p1, p2))
}



## LPI FUNCTIONS ## ------------------------------------------------------------

calclpi_exp <- function(Ni_init = N0i, Nj_init = N0j, n_pairedpops = 1, timesteps = steps, Ri = r_i, Rj = r_j, int_ij = alpha_ij, int_ji = alpha_ji, K){
  
  # initialize vector to store results (population sizes)
  Ni <- Ni_init
  Nj <- Nj_init
  
  # Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
  # calculate population sizes
  for(t in 2:timesteps){
    
    # population i
    temp_i = Ni[t-1]*(1 + Ri[1,t]*(1 - (Ni[t-1] + int_ij*Nj[t-1])/K[t]))
    temp_i[which(temp_i < 0)] <- 0 # assign 0 to negative population sizes
    Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
    
    # population j
    temp_j = Nj[t-1]*(1 + Rj[1,t]*(1 - (Nj[t-1] + int_ji*Ni[t-1])/K[t]))
    temp_j[which(temp_j < 0)] <- 0 # assign 0 to negative population sizes
    Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
  }
  
  # initialize df to store dts
  lambda_df = data.frame(dt_i = 0, dt_j = 0)
  for(i in 2:steps){
    # calculate dt
    lambda_df[i, "dt_i"] = log10(Ni[i]/Ni[i-1])
    lambda_df[i, "dt_j"] = log10(Nj[i]/Nj[i-1])
  }
  # take geometric mean
  lambda_gm <- apply(10^(lambda_df), 1, gm_mean)
  
  # calculate LPI
  lpi_exp <- data.frame(time = 1:steps, LPI = calclpi(log10(lambda_gm)))
  
  return(list("Ni" = Ni, "Nj" = Nj, "dt" = lambda_df, "dt_mean" = lambda_gm, "lpi" = lpi_exp))
}


# function to run GAM on population time series + predict over all time steps
gam_lpi <- function(pops){
  
  # estimate the smoothing parameter to be half the time series length
  smoothParm = round(length(unique(pops$time))/2)
  
  # transform population sizes
  pops$N <- log10(pops$N + 1)
  # wrangle into wide format
  pops_w = pivot_wider(data = pops, id_cols = time, names_from = popID, values_from = N) %>%
    mutate_at(vars(time), as.integer) %>% as.matrix()
  
  # run GAM (modified code from CalcLPI function in rlpi package)
  m <- list()
  for(i in 2:ncol(pops_w)){
    N <- as.vector(pops_w[,i])
    time <- as.vector(pops_w[,"time"])
    m[[i-1]] <- gam(N ~ s(time, k = smoothParm), 
                    family = gaussian(), fx = TRUE, method = "REML")
  }
  names(m) <- unique(pops$popID)
  
  # predict over time period
  pred_ls = lapply(m, predict.gam, type = "response", se.fit = TRUE)
  
  # wrangle into one long format dataframe
  pred = pred_ls %>%
    lapply(bind_cols) %>%
    lapply(mutate, time = pops_w[,"time"]) %>%
    bind_rows(.id = "popID")
  # join to observation dataframe
  pred = full_join(pops, pred, by = c("popID", "time"))
  
  return(pred)
}


# function to calculate dt from GAM predictions
get_dt <- function(gam.pred_ls, time){
  # gam.pred_ls: list of predictions from GAMs (one per population)
  
  # extract predicted population size values
  N = gam.pred_ls$fit 
  # assign errors to values for propagation
  errors(N) = gam.pred_ls$se.fit
  # un-log10
  N = 10^N # un-log
  
  # initialize df to store dts
  dt_df = data.frame(time = time, dt = NA, se = NA)
  for(i in 2:length(N)){
    # calculate dt
    dt = log10(N[i]/N[i-1])
    dt_df[i, "dt"] = dt # save in the table
    # save propagated error 
    dt_df[i, "se"] = unlist(errors(dt))
  }
  return("dt" = dt_df)
}

# geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# function to calculate geometric mean that will be used in bootstrapping
geoMean_boot <- function(dt, i){
  d <- dt[i] # indexing happens here
  avg <- gm_mean(d, na.rm = TRUE)
  return(avg)
}

# function for bootstrap confidence interval 
dt_boot = function(dt){
  # bootstrap resampling
  dt_r = boot::boot(dt, statistic = geoMean_boot, R = 1000)
  # calculate 95% confidence intervals
  dt_ci = boot::boot.ci(dt_r, type = "basic", R = 1000)
  # wrangle into table for output
  dtboot_df = data.frame(gm = dt_ci$t0, cilo = dt_ci$basic[4], cihi = dt_ci$basic[5])
  return(dtboot_df)
}

# function to calculate LPI from mean dt
calclpi <- function(dt){
  # calculate index value
  lpi = c(1) # initial value is 1 
  for(i in 2:length(dt)){
    lpi[i] <- lpi[i-1]*10^dt[i] }
  return(lpi)
}