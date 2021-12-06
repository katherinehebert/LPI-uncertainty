# Script to compute different measures of precision, including width of confidence intervals,
# success rates for capturing the ture trend within confidence intervals, etc.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

# function to compute these metrics for each scenario
compare_precision <- function(simID){
  
  # from chain calculation on raw abundances
  dt_chain <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/scenario", simID, "_l_dtchain.RDS"))
  dt_chain$time <- (dt_chain$time + 1969) %>% as.character()
  
  # from bootstrapped LPI values 
  bootstrap <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/rlpi/scenario", simID, "_bootstrap_rlpi.RDS"))
  
  # import true LPI values (no error)
  true <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/scenario", simID, "_true_rlpi.RDS"))
  true$time = as.character(true$time)
  
  # calculate LPI for each population
  
  # function to calculate LPI from mean dt
  calclpi <- function(dt){
    # calculate index value
    lpi = c(1) # initial value is 1 
    for(i in 2:length(dt)){
      lpi[i] <- lpi[i-1]*10^dt[i] }
    return(lpi)
  }
  
  # convert dt df to wide format
  dt_chain <- dt_chain %>% group_split(popID) %>%
    lapply(function(x) {
      x$lpi <- calclpi(x$dt_chain) 
      return(x)}
    ) %>% bind_rows()
  
  # convert to wide format 
  chain <- subset(dt_chain, select = c(time, lpi, popID)) %>%
    pivot_wider(names_from = time, values_from = lpi) %>%
    subset(select = -popID)
  # match up the time steps
  bootstrap <- bootstrap[,1:11]
  colnames(bootstrap) <- colnames(chain)
  
  # save the wide formats 
  bootstrap_w <- bootstrap
  chain_w <- chain
  
  ## plot distribution of LPI values ####
  
  # convert to long format for plotting
  bootstrap <- as.data.frame(bootstrap) %>% 
    pivot_longer(cols = everything(), names_to = "year", values_to = "bootstrap")
  chain <- as.data.frame(chain) %>% 
    pivot_longer(cols = everything(), names_to = "year", values_to = "chain")
  
  # join them together
  dt_compare_w <- inner_join(bootstrap, chain, by = "year") 
  dt_compare <- pivot_longer(dt_compare_w, cols = c(bootstrap, chain), names_to = "type", values_to = "lpi")
  
  # per year #### IMPORTANT ONE #######
  ggplot(dt_compare, aes(y = year, x = lpi, lty = type, fill = factor(stat(quantile)))) +
    stat_density_ridges(geom = "density_ridges_gradient", 
                        alpha = .7, lwd = 0, 
                        calc_ecdf = TRUE,
                        quantiles = c(0.025, 0.975),
                        scale = .9) +
    stat_density_ridges(aes(fill = NA), alpha = 0, lwd = .3, quantiles = c(.5), quantile_lines = TRUE,
                        scale = .9) +
    scale_fill_manual(
      name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
    ) + 
    scale_linetype_manual(values = c(1, 5),
                          name = "", 
                          labels = c("smoothed LPI", "raw LPI")
    ) + 
    labs(x = "LPI values", y = "Year", fill = "", title = paste0("Scenario ", simID))
  ggsave(paste0("figures/scenario", simID, "_precisiondensity.png"), width = 8.6, height = 7)
  
  
  # initialise a data frame to store results
  df <- data.frame(
    "time" = colnames(bootstrap_w),
    "mean_rlpi" = apply(bootstrap_w, 2, mean),
    "mean_chain" = apply(chain_w, 2, mean),
    "sd_rlpi" = apply(bootstrap_w, 2, sd),
    "sd_chain" = apply(chain_w, 2, mean),
    "q025_rlpi" = apply(bootstrap_w, 2, quantile, probs = .025),
    "q500_rlpi" = apply(bootstrap_w, 2, quantile, probs = .5),
    "q975_rlpi" = apply(bootstrap_w, 2, quantile, probs = .975),
    "q025_chain" = apply(chain_w, 2, quantile, probs = .025),
    "q500_chain" = apply(chain_w, 2, quantile, probs = .5),
    "q975_chain" = apply(chain_w, 2, quantile, probs = .975)
  )
  # join with true lpi dataset
  df <- left_join(df, subset(true, select = c(time, LPI_final)))
  df <- rename(df, "LPI_final_true" = "LPI_final")
  
  # do the CI capture the underlying mean/median/truth?
  df$mean_within_rlpiCI <- NA
  df$median_within_rlpiCI <- NA
  for(i in 1:nrow(df)){
    # determine whether chain lpi value is within the LPI intervals
    if(df$mean_chain[i] >= df$q025_rlpi[i] & df$mean_chain[i] <= df$q975_rlpi[i]){  
      df$mean_within_rlpiCI[i] = "Success"
    } else df$mean_within_rlpiCI[i] = "Failure"
    # repeat but with median
    if(df$q500_chain[i] >= df$q025_rlpi[i] & df$q500_chain[i] <= df$q975_rlpi[i]){  
      df$median_within_rlpiCI[i] = "Success"
    } else df$median_within_rlpiCI[i] = "Failure"
  }
  
  # determine whether the true value is within the LPI
  df$true_within_rlpiCI <- NA
  for(i in 1:nrow(df)){
    # determine whether chain lpi value is within the LPI intervals
    if(df$LPI_final_true[i] >= df$q025_rlpi[i] & df$LPI_final_true[i] <= df$q975_rlpi[i]){  
      df$true_within_rlpiCI[i] = "Success"
    } else df$true_within_rlpiCI[i] = "Failure"
  }
  
  # determine difference between mean and the CI limits
  df$mean_diff_rlpiCI_low <- NA
  df$mean_diff_rlpiCI_high <- NA
  for(i in 1:nrow(df)){
    if(df$mean_within_rlpiCI[i] == "Failure"){
      if(df$mean_chain[i] < df$q025_rlpi[i]){
        df$mean_diff_rlpiCI_low[i] <- df$q025_rlpi[i] - df$mean_chain[i]
      } else if(df$mean_chain[i] > df$q975_rlpi[i]) {
        df$mean_diff_rlpiCI_high[i] <- df$mean_chain[i] - df$q975_rlpi[i] 
      }
    } else next
  }
  # if yes, yes
  # if no, is it below the CI_low? if yes, by how much?
  # if no, is it above the CI_high? if yes, by how much?
  # determine difference between mean and the CI limits
  df$true_diff_rlpiCI <- 0
  for(i in 1:nrow(df)){
    if(df$true_within_rlpiCI[i] == "Failure"){
      if(df$LPI_final_true[i] < df$q025_rlpi[i]){
        df$true_diff_rlpiCI[i] <- df$LPI_final_true[i] - df$q025_rlpi[i]  
      } else if(df$LPI_final_true[i] > df$q975_rlpi[i]) {
        df$true_diff_rlpiCI[i] <- df$LPI_final_true[i] - df$q975_rlpi[i] 
      }
    } else next
  }
  
  # calculate interval widths & difference in interval width
  df$rlpiCI_width <- df$q975_rlpi - df$q025_rlpi
  df$chainCI_width <- df$q975_chain - df$q025_chain
  df$rlpiCI_relativewidth <- (df$q975_rlpi - df$q025_rlpi)/df$LPI_final
  df$chainCI_relativewidth <- (df$q975_chain - df$q025_chain)/df$mean_chain
  df$CI_diff <- (df$rlpiCI_width - df$chainCI_width)/df$LPI_final 
  
  
  # calculate overlap between distribution of chain dt and rlpi dt
  df$overlap = group_split(dt_compare_w, year) %>%
    lapply(
      function(x){
      if(x$year[1] == "1970"){
        return(NA)
        } else { 
        bayestestR::overlap(x$bootstrap, x$chain)[1] 
  }
        }) %>% unlist()
  
  
  # calculate percentile of the true mean dt within the bootstrapped dt distribution
  temporary <- group_split(dt_compare_w, year) 
  df$percentile = NA
  for(i in 2:length(temporary)){ # skip first time step bc baseline
    percentile <- ecdf(temporary[[i]]$bootstrap)
    df$percentile[i] <- percentile(true$LPI_final[i])
  }
  
  saveRDS(df, paste0("outputs/scenario", simID, "_precision.RDS"))
  
}

# get all scenario names
sim_ids <- c(
  paste0("1", LETTERS[1:9]),
  paste0("2", LETTERS[1:18]),
  paste0("3", LETTERS[1:18]),
  paste0("4", LETTERS[1:18]),
  paste0("5", LETTERS[1:18]),
  paste0("6", LETTERS[1:18]),
  paste0("7", LETTERS[1:18])
)
# run for all scenarios
lapply(sim_ids, compare_precision)


# #### CALC SUCCESS RATE - this should go in accuracy_and_precision.Rmd.
# #### Maybe plot across multiple scenarios rather than at each step?
# df %>% 
#   group_by(direction) %>% 
#   summarise(p_success = length(which(true_within_rlpiCI == "Success"))/length(true_within_rlpiCI))
# ##############

# do the GAM errors correspond to the rror introduced into the time series?

# import true LPI values (no error)
# get all scenario names
sim_ids <- c(
  paste0("1", LETTERS[1:9]),
  paste0("2", LETTERS[1:18]),
  paste0("3", LETTERS[1:18]),
  paste0("4", LETTERS[1:18]),
  paste0("5", LETTERS[1:18]),
  paste0("6", LETTERS[1:18]),
  paste0("7", LETTERS[1:18])
)
gams <- lapply(sim_ids, function(x) readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/models/rlpi/scenario", x, ".rds")))
se <- lapply(gams, predict, se.fit = TRUE)

gamse <- list()
for(i in 1:length(se)) {
  gamse[[i]] <- data.frame("se" = se[[i]]$se.fit)
  colnames(gamse[[i]]) <- sim_ids[i]
}
gamse <- bind_cols(gamse)

pdf("outputs/measurementerror_vs_GAMerror.pdf")

boxplot(gamse[,grep("1", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen")) # attempt
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("2", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen"))
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("3", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen"))
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("4", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen"))
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("5", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen"))
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("6", colnames(gamse))], ylim = c(0, 0.1), 
        col = c("tomato2", "dodgerblue2", "mediumseagreen"))
abline(h = 0.05, lty = 2)

boxplot(gamse[,grep("7", colnames(gamse))], ylim = c(0, 0.1))
abline(h = 0.05, lty = 2)

dev.off()