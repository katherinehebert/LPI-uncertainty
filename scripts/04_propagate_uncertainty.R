# Script to propagate uncertainty in the calculation of the LPI
# following equations 1-8 by Dominique Gravel

library(dplyr)
library(tidyr)

# get scenario IDs
scenarios <- gsub("_l.RDS", "", list.files(path = "simulations/", pattern = "_l.RDS")) %>%
  gsub("scenario", "", .)

# calculate uncertainty in each scenario

for(scenarioID in scenarios){

  ## load simulated populations --------------------------------------------------------
  
  Nraw <- readRDS(paste0("simulations/scenario", scenarioID, "_l.RDS"))
  # add 1 to any extinctions to avoid dividing by 0
  if(!is.null(length(which(Nraw$N == 0)))){Nraw$N[which(Nraw$N == 0)] = 1}
  Nparams <- readRDS(paste0("simulations/scenario", scenarioID, "_params.RDS"))
  
  # get sigmas
  sigma_m <- dplyr::filter(Nparams, parameter == "Observation error")[,2:3] %>% as.vector() %>% apply(1:2, function(x) x^2) %>% as.vector()
  sigma_p <- dplyr::filter(Nparams, parameter == "Process error")[,2:3] %>% as.vector() %>% apply(1:2, function(x) x^2) %>% as.vector()
  
  # pivot to wide format 
  N <- subset(Nraw, select = c(popID, time, N)) %>% 
    tidyr::pivot_wider(names_from = time, 
                       values_from = N)
  rownames(N) <- N$popID
  N <- subset(N, select = -c(popID))
  
  #### estimate population growth rate trend #### --------------------------------
  
  # equation 2 ----
  # expectation of the growth rate trend
  eq2 <- function(N, sigma_measure = 0){
    d <- c()
    for(t in 2:length(N)){
      d[t] <- log10(N[t]/N[t-1]) + (sigma_measure^2)/(2*(N[t-1]^2 - N[t]^2))
    }
    return(d)
  }
  
  dt <- cbind(apply(N[1:10,], 1, eq2, sigma_measure = sigma_m[1]),
              apply(N[11:20,], 1, eq2, sigma_measure = sigma_m[2]))
  
  
  # equation 2
  # just the growth rate part
  eq2_donly <- function(N){
    eq2dt <- c()
    for(t in 2:length(N)){
      eq2dt[t] <- log10(N[t]/N[t-1])
    }
    return(eq2dt)
  }
  donly <- apply(N, 1, eq2_donly)
  
  
  # equation 2
  # just the uncertainty part
  eq2_varonly <- function(N, sigma_measure = 0){
    eq2var <- c()
    for(t in 2:length(N)){
      eq2var[t] <- (sigma_measure^2)/(2*((N[t-1])^2 - (N[t])^2))
      # should the population sizes be logged here?
      #eq2var[t] <- (sigma_measure^2)/(2*(log(N[t-1])^2 - log(N[t])^2))
    }
    return(eq2var)
  }
  varonly <- cbind(apply(N[1:10,], 1, eq2_varonly, sigma_measure = sigma_m[1]),
                  apply(N[11:20,], 1, eq2_varonly, sigma_measure = sigma_m[2]))
  
  save_this <- list("dt" = dt, "d" = donly, "measerr_correction" = varonly)
  saveRDS(save_this, paste0("outputs/poptrend_estimate/scenario", scenarioID, ".rds"))
  
  # par(mfrow = c(2,2))
  # matplot(t(N), type = "l", main = "Population size", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # matplot(dt, type = "l", main = "Growth rate (eq2)", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # matplot(donly, type = "l", main = "Growth rate without uncertainty", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # matplot(varonly, type = "l", main = "Measurement uncertainty", add = TRUE, col = PNWColors::pnw_palette("Sunset2", 5), lty= 2)
  # matplot(varonly, type = "l", main = "Measurement uncertainty", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  
  ## equation 3 ----
  
  eq3 <- function(N, sigma_measure = 0, sigma_process = 0){
    eq3var <- c()
    for(t in 2:length(N)){
      eq3var[t] <- sigma_process^2 + (sigma_measure^2)*((1/((N[t])^2) - (1/(N[t-1])^2)))
    }
    return(eq3var)
  }
  var_dt <- cbind(apply(N[1:10,], 1, eq3, sigma_measure = sigma_m[1], sigma_process = sigma_p[1]),
                   apply(N[11:20,], 1, eq3, sigma_measure = sigma_m[2], sigma_process = sigma_p[2]))
  
  # par(mfrow = c(1,3))
  # matplot(t(N), type = "l", main = "Population size", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # matplot(dt, type = "l", main = "Growth rate (eq2)", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # matplot(var_dt, type = "l", main = "Uncertainty in growth rate (eq3)", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  
  
  #### average growth rate #### --------------------------------------------------
  
  ## equation 4 ----
  ## average growth rate at each time step
  
  dt_bar <- apply(dt, 1, mean)
  
  # par(mfrow = c(1,1))
  # matplot(dt, type = "l", main = "Growth rate (eq3)", col = PNWColors::pnw_palette("Sunset2", 5), lty= 1)
  # lines(dt_bar, lwd = 3)
  
  ## equation 5 ----
  ## uncertainty in the average growth rate
  
  dt_cov <- cov(dt, use = "pairwise.complete.obs")
  dt_cov <- dt_cov[which(lower.tri(dt_cov))]
  
  var_dtbar = (1/nrow(N))*(apply(var_dt, 1, sum) + 2*(abs(sum(dt_cov)))) # FLAG ---- 
  
  #plot(var_dtbar, type = "l")
  
  ## equation 6 ----
  ## calculate LPI (without uncertainty correction)
  
  # function to calculate LPI value without uncertainty correction
  calclpi <- function(dt_bar){
    I = 1 
    for(i in 2:length(dt_bar)){
      I[i] <- I[i-1]*10^dt_bar[i]
    }
    return(I)
  }
  
  ## equation 7 ----
  
  # function to calculate LPI value WITH uncertainty correction
  calclpi_corrected <- function(dt_bar){
    I = 1 
    for(i in 2:length(dt_bar)){
      I[i] <- I[i-1]*10^dt_bar[i] + 0.5*(10^dt_bar[i]*var_dtbar[i])
    }
    return(I)
  }
  
  lpi_nocorrection = calclpi(dt_bar)
  lpi_correction = calclpi_corrected(dt_bar)
  
  plot(lpi_nocorrection, type = "l", ylab = "I")
  lines(lpi_correction, col = "purple")
  
  ## equation 8 ----
  
  # function to obtain the variance of the LPI
  
  var_lpi <- (10^(2*dt_bar))*var_dtbar
  #plot(var_lpi, type = "l")
  
  
  # results to save
  
  lpi_res <- data.frame(
    "time" = 1:11,
    "dtbar" = dt_bar,
    "dtbar_variance" = var_dtbar,
    "lpi_nocorrection" = lpi_nocorrection, # eq 6
    "lpi_correction" = lpi_correction, # eq 7
    "lpi_variance" = var_lpi, # eq 8,
    "lpi_bias" = lpi_nocorrection - lpi_correction
  )
  saveRDS(lpi_res, paste0("outputs/scenario", scenarioID, "_uncertaintypropagation.RDS"))
}


## check out the results ------

## load the results
res <- lapply(paste0("outputs/", list.files("outputs/", "_uncertaintypropagation.RDS")[-c(1,2)]), readRDS)
names(res) <- scenarios
res <- bind_rows(res, .id = "scenario")
res$scenario <- paste0("scenario", res$scenario)
saveRDS(res, "outputs/all_uncertaintypropagation.RDS")
