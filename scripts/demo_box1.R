# Script to demonstrate the need for a multispecies index of biodiversity change
# using a simple ecological example

library(tidyverse)
library(hrbrthemes)
library(patchwork)
library(synchrony)
library(errors)

# load required functions
source('scripts/scenario_functions.R')

# 1. explore data to find a good example ----

# read LPD dataset
#lpd <- read_csv("data_raw/LPR2020data_public.csv")

# subset to example system
ex_sub <- lpd %>% 
  filter(Location == "Kluane Lake Research Station, Yukon, Canada") 

# clean up and convert to long format
ex <- ex_sub %>%
  pivot_longer(cols = "1950":ncol(.), values_to = "size", names_to = "year") %>%
  mutate_at(vars(year), as.integer) %>%
  mutate_at(vars(size), as.numeric) %>%
  # remove rows with NA
  drop_na(size)

# create clean wide version
ex_wide <- ex_sub[,30:ncol(ex_sub)] %>% 
  apply(1:2, as.numeric) %>%
  t()
colnames(ex_wide) = ex_sub$Binomial

# plot time series ----
p_ts <- ggplot(ex) +
  geom_line(aes(x = year, y = size, col = Binomial), lwd = 1.1) +
  scale_y_log10() +  
  coord_cartesian(xlim = c(1990, max(ex$year))) + # crop 
  theme_ipsum_rc() +
  labs(y = "Population size (log10)"
       #title = "Interacting populations",
       #subtitle = "Coyote, Snowshoe hare, and Canadian lynx"
       )


# calculate synchrony ----

trio_sync <- community.sync(na.omit(ex_wide), nrands = 100)
plot(trio_sync)

ex_sync <- as.data.frame(ex_wide)
ex_sync$year <- rownames(ex_wide)
ex_sync <- relocate(.data = ex_sync, year) %>% 
  mutate_at(vars(year), as.integer) %>%
  as.matrix() %>% na.omit()

sync_canlep <- phase.sync(t1 = ex_sync[,c(1,2)], 
                          t2 = ex_sync[,c(1,3)],
                          nrands = 100)
plot(sync_canlep$deltaphase$phasediff ~ sync_canlep$deltaphase$timestep, pch = 16)
# redo this ^^ but with growth rates....

# get step-wise growth rates ----

# function to calculate population growth rate (dt) using the chain method
# N = vector of population sizes (not transformed)
calc_dt_chain <- function(N){
  dt = c(0) # initial value
  for(i in 2:length(N)){
    dt[i] = log10(N[i]/N[i-1])
  }
  return(dt)
}

# calculate growth rate
dt <- ex_sync
for(i in 2:4){
  dt[,i] <- calc_dt_chain(dt[,i])
}
ex_dt <- as.data.frame(dt) %>%
  pivot_longer(cols = 2:4, names_to = "Binomial", values_to = "size")

# plot growth rates ----
p_dt <- ggplot(ex_dt) +
  geom_line(aes(x = year, y = size, col = Binomial), lwd = 1.1) +
  #scale_y_log10() +  
  theme_ipsum_rc() +
  labs(y = "Growth rate (log10)"
       #title = "Interacting populations",
       #subtitle = "Coyote, Snowshoe hare, and Canadian lynx"
  )

(p_ts / p_dt)

# calculate linear covariance between populations ----

cov_ex_dt <- dt[,2:4] %>% cov() %>% cov2cor()
corrplot::corrplot(cov_ex_dt, type = "upper", method = "color",
                   addCoef.col = "white", # Add coefficient of correlation
                   tl.col="black", tl.srt = 45)

# calculate synchrony ----

trio_sync <- community.sync(dt[,2:4], nrands = 100)
plot(trio_sync)


#### mean growth rate per time step - no cov

mean_dt <- data.frame(
  year = dt[,1],
  dt = apply(10^dt[,2:4], 1, EnvStats::geoMean),
  sd = apply(10^dt[,2:4], 1, EnvStats::geoSD)
)

# plot mean growth rate ----
p_mean <- ggplot(mean_dt) +
  geom_ribbon(aes(x = year, 
                  ymin = dt - 1.96*sd,
                  ymax = dt + 1.96*sd)) +
  geom_line(aes(x = year, y = dt), color = "white", lwd = 1.1) +
  #scale_y_log10() +  
  theme_ipsum_rc() +
  labs(y = "Growth rate (log10)",
       title = "Mean growth rate (without covariance)")
p_mean


#### GAM ----

# run GAM (modified code from CalcLPI function in rlpi package)
m <- list()
for(i in 2:ncol(dt)){
  N <- as.vector(dt[,i])
  time <- as.vector(dt[,"year"])
  m[[i-1]] <- mgcv::gam(N ~ s(time, k = round(nrow(dt)/2)), 
                  family = gaussian(), fx = TRUE, method = "REML")
}

# predict over time period
names(m) <- colnames(dt)[2:4]
gams <- lapply(m, mgcv::predict.gam, type = "response", se.fit = TRUE) %>%
  bind_rows(.id = "Binomial") %>%
  mutate(year = rep(dt[,1], 3)) %>%
  as.data.frame()

# plot mean growth rate ----
p_gams <- ggplot(gams) +
  geom_ribbon(aes(x = year,
                  ymin = fit - 1.96*se.fit,
                  ymax = fit + 1.96*se.fit)) +
  geom_line(aes(x = year, 
                y = fit), 
            color = "white", lwd = 1.1) +
  #scale_y_log10() +  
  theme_ipsum_rc() +
  facet_wrap(~Binomial) +
  labs(y = "Growth rate (log10)")
p_gams

