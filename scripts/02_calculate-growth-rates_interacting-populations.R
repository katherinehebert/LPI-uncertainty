# Script to calculate growth rates of simulated populations from GAMs (as in LPI)

# load required packages
require(mgcv)
require(tidyverse)
require(errors)
# set ggplot theme
theme_set(theme_linedraw() + theme(panel.grid = element_blank()))


# import simulated populations
pops <- readRDS("simulations/paired_antagonistic_l.RDS") %>%
  mutate_at(vars(time), as.integer)

# 1. GAMs ----------------------------------------------------------------------

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

# plot the predicted trend over the original
ggplot(pred, aes(x = time, group = popID)) +
  geom_line(aes(y = N, col = pop), lty = 2, lwd = .5) +
  geom_ribbon(aes(ymin = fit - se.fit,
                  ymax = fit + se.fit, fill = pop), alpha = .3) +
  geom_line(aes(y = fit)) +
  labs(y = "N(t) (log10)", col = "Paired \npopulations", fill = "Paired \npopulations",
       caption = "Ribbon shows standard error from GAM predictions.") + 
  facet_wrap(~set) 
ggsave(filename = "paired_antagonistic_gampredictions.png", path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")

# 2. Growth rates from GAMs ----------------------------------------------------

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

# calculate dt from GAM predictions
dt <- lapply(pred_ls, get_dt, time = pops_w[,"time"]) %>% bind_rows(.id = "popID")
# join to observations dataframe
dt_df = full_join(pops, dt, by = c("popID", "time"))

# plot dts
ggplot(data = dt_df, aes(x = time, group = popID)) +
  geom_line(aes(y=dt, col = pop)) + 
  geom_ribbon(aes(ymin = dt-se, ymax = dt+se, fill = pop), alpha = .2) +
  facet_wrap(~set) + 
  labs(y = "Growth rate (log10)", 
      col = "Paired \npopulations", fill = "Paired \npopulations",
      caption = "Ribbon shows propagated standard error from GAM predictions.") +
  theme(legend.position = "none")

# save outputs -----------------------------------------------------------------
saveRDS(m, "outputs/paired_antagonistic_GAM.RDS")
saveRDS(dt_df, "outputs/paired_antagonistic_dt_l.RDS")
ggsave(filename = "paired_antagonistic_dt.png", path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")
