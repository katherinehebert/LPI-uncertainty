# Function to calculate the expected and observed Living Planet Index from 
# growth rates of paired interacting populations and measuring uncertainty via 
# bootstrapping and error propagation.

# simname = name used for the outputs and figures for this simulation
# setID = identifying character for the set of interacting populations (e.g. i and j)

get_lpi <- function(simname){

# import growth rates
dt_df <- readRDS(paste0("outputs/", simname, "_results.RDS"))

# calculate expected lpi
truth <- readRDS(paste0("outputs/", simname, "_true.RDS"))
temp <- truth %>% group_by(time) %>% group_split()
temp <- lapply(temp, function(x) gm_mean(10^x$dt)) %>% unlist() %>% unname() %>% log10()
lpi_exp <- calclpi(dt = temp)

# 1. Bootstrapping -------------------------------------------------------------

# calculate geometric mean growth rate with 95% CI from bootstrapping
dt_gm <- list(data.frame(gm = 1, cilo = 1, cihi = 1))
for(t in 2:length(unique(dt_df$time))){
  dt_v = dt_df[which(dt_df$time == t), "dt"] %>% na.omit()
  dt_gm[[t]] = dt_boot(10^dt_v$dt) 
}
dt_gm = bind_rows(dt_gm) %>% log10() %>% 
  mutate(time = unique(dt_df$time))


# Calculate LPI with bootstrapped 95%CI  -----------

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
for(i in 2:length(dt_ls)){
  # extract growth rates and assign error
  dt_err = dt_ls[[i]]$dt
  errors(dt_err) = dt_ls[[i]]$se 
  # calculate geometric mean
  gm <- gm_mean(dt_err) %>% log10()
  # wrangle into a vector (mean, error)
  dt_gm_err[i, "gm"] = drop_errors(gm)
  dt_gm_err[i, "err"] = errors(gm)
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
lpi$LPI_true <- lpi_exp
saveRDS(lpi, paste0("outputs/", simname, "_lpi.RDS"))

# Plot LPI with both errors  ------------
ggplot(lpi, aes(x = time)) +
  # bootstrap version
  geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot), alpha = .2, fill = "blue") +
  geom_line(aes(y = LPI_boot)) +
  # propagation version
  geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se), alpha = .2, fill = "#e7298a") +
  geom_line(aes(y = LPI_se)) +
  # expected 
  ylim(c(0,2)) + labs(y = "Living Planet Index", x = "")
ggsave(filename = paste0(simname, "_lpi.png"), path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")

return(lpi)
}