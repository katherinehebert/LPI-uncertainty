# Script to calculate Living Planet Index of simulated populations 
# with error intervals calculated via bootstrapping

# load required packages
require(tidyverse)
require(errors)
# set ggplot theme
theme_set(theme_linedraw() + theme(panel.grid = element_blank()))


# import growth rates
dt_df <- readRDS("outputs/paired_antagonistic_dt_l.RDS")

# geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# 1. Bootstrapping -------------------------------------------------------------

  # Get geometric mean growth rate ----------------
  
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
  
  # calculate geometric mean growth rate with 95% CI from bootstrapping
  dt_gm <- list(data.frame(gm = 1, cilo = 1, cihi = 1))
  for(t in 2:length(unique(dt_df$time))){
    dt_v = dt_df[which(dt_df$time == t), "dt"] %>% na.omit()
    dt_gm[[t]] = dt_boot(10^dt_v$dt) 
  }
  dt_gm = bind_rows(dt_gm) %>% log10() %>% 
    mutate(time = unique(dt_df$time))
  
  
  # Calculate LPI with bootstrapped 95%CI  -----------
  
  # function to calculate LPI from mean dt
  calclpi <- function(dt){
    # calculate index value
    lpi = c(1) # initial value is 1 
    for(i in 2:length(dt)){
      lpi[i] <- lpi[i-1]*10^dt[i] }
    return(lpi)
  }
  
  # calculate LPI with confidence intervals
  lpi_boot <- data.frame(
    time = 1:nrow(dt_gm), 
    LPI_boot = calclpi(dt_gm$gm),
    cilo_boot = calclpi(dt_gm$cilo),
    cihi_boot = calclpi(dt_gm$cihi)
  )
  

# 2. Error propagation ---------------------------------------------------------

# prepare growth rates
dt_ls <- mutate_at(dt_df, vars(dt), function(x) 10^x) %>%
  mutate_at(vars(time), as.factor) %>%
  group_by(time) %>% group_split()
# remove first time step (all NA because there was no previous value) 
dt_ls[[1]] <- NULL

# take geometric mean growth rate per time step (with error propagation)
dt_gm_err <- data.frame("time" = 1, "gm" = 1, "err" = 0)
for(i in 1:length(dt_ls)){
  # extract growth rates and assign error
  dt_err = dt_ls[[i]]$dt
  errors(dt_err) = dt_ls[[i]]$se 
  # calculate geometric mean
  gm <- gm_mean(dt_err) %>% log10()
  # wrangle into a vector (mean, error)
  dt_gm_err[i+1, "gm"] = drop_errors(gm)
  dt_gm_err[i+1, "err"] = errors(gm)
  }  
dt_gm_err$time <- 1:nrow(dt_gm_err)


# calculate LPI with %95CI from propagated error
lpi_err <- data.frame(
  time = dt_gm_err$time, 
  LPI_se = calclpi(dt_gm_err$gm),
  cilo_se = calclpi(dt_gm_err$gm - 1.96*dt_gm_err$err),
  cihi_se = calclpi(dt_gm_err$gm + 1.96*dt_gm_err$err)
)

# 3. Plot to compare approaches ------------------------------------------------

# join as one df
lpi <- inner_join(lpi_boot, lpi_err, by = "time")

# Plot LPI with bootstrapped 95%CI  ------------
ggplot(lpi, aes(x = time)) +
  # bootstrap version
  geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot), alpha = .2, fill = "blue") +
  geom_line(aes(y = LPI_boot)) +
  # propagation version
  geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se), alpha = .4, fill = "#e7298a") +
  geom_line(aes(y = LPI_se)) +
  # expected 
  geom_hline(aes(yintercept = 1), lty = 2) +
  # annotations
  annotate("text", x = 9, y = 1.3, label = "Bootstrap", 
           colour = "darkblue", fontface = "italic", angle = 8) +
  annotate("text", x = 9, y = 0.75, label = "Propagation", 
           colour = "#980043", fontface = "italic", angle = -15) +
  ylim(c(0.5,1.5)) + labs(y = "Living Planet Index")

# save outputs -----------------------------------------------------------------
saveRDS(lpi, "outputs/paired_antagonistic_LPI.RDS")
ggsave(filename = "paired_antagonistic_LPI.png", path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")
