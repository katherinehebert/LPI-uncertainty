# Script to gather all results into one table for plotting

# load packages and prep environment
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
theme_set(ggpubr::theme_pubr())


# loading outputs ----

## MY LPI ##

# import results of each scenario's LPI
lpi <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_lpi.RDS")), readRDS)
names(lpi) <- gsub("_lpi.RDS", "", list.files(path = "outputs/", pattern = "_lpi.RDS"))
# and bind into one data frame
lpi <- dplyr::bind_rows(lpi, .id = "scenario")

## RLPI ##

# import results from rlpi package
temp = list.files(path = "outputs/", pattern = "_rlpi.RDS")[-grep("true", list.files(path = "outputs/", pattern = "_rlpi.RDS"))]
rLpi <- lapply(paste0("outputs/", temp), readRDS)
names(rLpi) <- gsub("_rlpi.RDS", "", temp)
# and bind into one data frame
rLpi <- dplyr::bind_rows(rLpi, .id = "scenario")
rLpi$time <- rLpi$time - 1969

# import truth results from rlpi package
temp = list.files(path = "outputs/", pattern = "_rlpi.RDS")[grep("true", list.files(path = "outputs/", pattern = "_rlpi.RDS"))]
rLpi_true <- lapply(paste0("outputs/", temp), readRDS)
names(rLpi_true) <- gsub("_true_rlpi.RDS", "", temp)
# and bind into one data frame
rLpi_true <- dplyr::bind_rows(rLpi_true, .id = "scenario")
rLpi_true$time <- rLpi_true$time - 1969
colnames(rLpi_true)[2:4] <- paste0(colnames(rLpi_true)[2:4], "_true")

# join the rlpi results together
rLpi = left_join(rLpi, rLpi_true)

## PARAMS ##

# import parameter tables for each scenario
params <- lapply(paste0("simulations/", list.files(path = "simulations/", pattern = "_params.RDS")), readRDS)
params <- lapply(params, data.table::transpose) %>%
  lapply(function(x){
    colnames(x) <- gsub(" ", "_", x[1,])
    x <- x[-c(1,2),]
  })
names(params) <- gsub("_params.RDS", "", list.files(path = "simulations/", pattern = "_params.RDS"))
# and bind into one data frame
params <- dplyr::bind_rows(params, .id = "scenario")

# import expected uncertainties of each scenario's growth rates
uncertainties <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_uncertainty.RDS")), readRDS)
names(uncertainties) <- gsub("_uncertainty.RDS", "", list.files(path = "outputs/", pattern = "_uncertainty.RDS"))
uncertainty <- bind_rows(uncertainties, .id = "scenario")

# make table for categories from carrying capacity scenarios
K_scenarios <- data.frame("scenario" = params$scenario, "direction" = NA)
K_scenarios$direction[which(K_scenarios$scenario %like% 'A|D|G|J|M|P')] <- "decline"
K_scenarios$direction[which(K_scenarios$scenario %like% 'B|E|H|K|N|Q')] <- "stable"
K_scenarios$direction[which(K_scenarios$scenario %like% 'C|F|I|L|O|R')] <- "growth"

# associate lag with the scenarios
params$Lag[which(params$scenario %like% 'scenario4|scenario5')] <- "1"
params$Lag[which(params$scenario %like% 'scenario6|scenario7')] <- "2"

# join all tables together
df <- dplyr::left_join(lpi, K_scenarios) %>% dplyr::left_join(params)
df <- left_join(df, rLpi)

# accuracy ----

# calculate LPI accuracy as % difference [(estimated - true)/true * 100]
#df$accuracy_boot <- ((df$LPI_boot - df$LPI_true)/df$LPI_true)*100
#df$accuracy_boot <- df$LPI_boot - df$LPI_true
df$accuracy_boot <- df$LPI_final - df$LPI_true

# precision ----

# does the true LPI fall within the 95% confidence interval?
if(df$LPI_true >= df$CI_low && df$LPI_true <= df$CI_high){
  df$precision_boot <- "yes"} else {
    df$precision_boot <- "no"
  }
# this next ifelse is untouched since correcting the LPI results with rlpi ---------
if(df$LPI_true <= df$cihi_se && df$LPI_true >= df$cilo_se){
  df$precision_se <- "yes"} else {
    df$precision_se <- "no"
  }

# interval width 
df$interval_width <- (df$CI_high-df$CI_low)/df$LPI_final #((df$cihi_boot-df$cilo_boot)/df$LPI_boot)
df$interval_width_se <- (df$cihi_se-df$cilo_se)/df$LPI_se #((df$cihi_se-df$cilo_se)/df$LPI_se)

# this next ifelse is untouched since correcting the LPI results with rlpi --------
# is the bootstrap interval within the error propagated one?
df$cilo_selargest <- df$cilo_se < df$cilo_boot
df$cihi_selargest <- df$cihi_se > df$cihi_boot
# put NA on the first time step, which has no uncertainty
df$cilo_selargest[which(df$time == 1)] <- NA
df$cihi_selargest[which(df$time == 1)] <- NA

# precision vs. uncertainty ----

# group uncertainties by scenario x time to match up with df
# temp <- uncertainty %>% group_split(scenario, time) 
# mean_uncertainty <- list()
# sd_uncertainty <- list()
# for(i in 1:length(temp)){
#   mean_uncertainty[[i]] <- mean(temp[[i]]$uncertainty_LPI)
#   sd_uncertainty[[i]] <- sd(temp[[i]]$uncertainty_LPI)
# }
# mean_uncertainty <- cbind(unlist(mean_uncertainty), unlist(sd_uncertainty))
# colnames(mean_uncertainty) <- c("mean_uncertainty", "sd_uncertainty")
# # attach to df
# df <- cbind(df, mean_uncertainty)

# calculate difference between the expected interval width and the bootstrapped one
#df$interval_diff <- (df$interval_width-df$mean_uncertainty)/df$LPI_true
#df$interval_diff <- (df$interval_width-df$mean_uncertainty)#/df$LPI_true
df$interval_diff <- (df$interval_width - df$interval_width_se)
# NOTE: what will this be now???? ---------------

# format df columns for plotting
df$Lag <- factor(df$Lag, levels = c("0", "1", "2"))
df$direction <- factor(df$direction, levels = c("decline", "stable", "growth"))
colnames(df)[11] <- "N0"
colnames(df)[12] <- "lambda"
colnames(df)[13] <- "interaction"

# save to file
saveRDS(df, "outputs/all_results.RDS")
readr::write_csv(df, "outputs/all_results.csv")