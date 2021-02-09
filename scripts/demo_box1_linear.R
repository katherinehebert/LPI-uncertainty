# Script to demonstrate the need for a multispecies index of biodiversity change
# using a simple ecological example

library(tidyverse)
library(hrbrthemes)
library(patchwork)
library(synchrony)
library(errors)
library(broom)

# load required functions
source('scripts/scenario_functions.R')

# set ggplot theme
theme_set(theme_ipsum_rc())

# 1. explore data to find a good example ----

# read LPD dataset
lpd <- read_csv("data_raw/LPR2020data_public.csv")

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

# reorganise data for the function
ex_sync <- as.data.frame(ex_wide)
ex_sync$year <- rownames(ex_wide)
ex_sync <- relocate(.data = ex_sync, year) %>%
  mutate_at(vars(year), as.integer) %>%
  as.matrix() %>% na.omit()

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

# plot dt(t) ~ d(t+1) - stability plot

dt_stability <- ex_dt %>%
  group_by(Binomial) %>%
  summarise(
    year = year,
    dt_0 = size,
    dt_1 = size,
    )

p_all <- ggplot(data = dt_stability) +
  geom_point(aes(x = dt_0, y = lag(dt_1), col = Binomial)) +
  geom_abline(aes(slope = 1, intercept = 0), lwd = .2)

p_each <- p_all + facet_wrap(~Binomial)

(p_all / p_each)

# calculate linear covariance between populations ----

cov_ex_dt <- dt[,2:4] %>% cov() %>% cov2cor()
corrplot::corrplot(cov_ex_dt, type = "upper", method = "color",
                   addCoef.col = "white", # Add coefficient of correlation
                   tl.col="black", tl.srt = 45)

# calculate synchrony ----

trio_sync <- community.sync(dt[,2:4], nrands = 100)
plot(trio_sync)


#### mean growth rate per time step (not treating populations as covarying)

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
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  #scale_y_log10() +  
  theme_ipsum_rc() +
  labs(y = "Growth rate (log10)")
p_mean


# simple linear model per time series

m_canis <- lm(dt[,2] ~ dt[,1])
pred_canis <- predict(m_canis, se = TRUE)
m_canis_summ <- tidy(m_canis)

m_lepus <- lm(dt[,3] ~ dt[,1])
pred_lepus <- predict(m_lepus, se = TRUE)
m_lepus_summ <- tidy(m_lepus)

m_lynx <- lm(dt[,4] ~ dt[,1])
pred_lynx <- predict(m_lynx, se = TRUE)
m_lynx_summ <- tidy(m_lynx)

# slopes
slopes_lm <- data.frame(
  population = c("Canis", "Lepus", "Lynx"),
  est = c(m_canis_summ$estimate[2],
          m_lepus_summ$estimate[2],
          m_lynx_summ$estimate[2]),
  se = c(m_canis_summ$std.error[2],
         m_lepus_summ$std.error[2],
         m_lynx_summ$std.error[2])
  )


preds_lm <- data.frame(
  population = rep(colnames(dt)[2:4], each = nrow(dt)),
  year = rep(dt[,1], 3),
  dt_lm = c(pred_canis$fit, pred_lepus$fit, pred_lynx$fit),
  se = c(pred_canis$se.fit, pred_lepus$se.fit, pred_lynx$se.fit)
)

ggplot(preds_lm, aes(x = year)) +
  geom_ribbon(aes(ymin = dt_lm - 1.96*se,
                  ymax = dt_lm + 1.96*se)) +
  geom_line(aes(y = dt_lm), col = "white") +
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  facet_wrap(~population) +
  labs(x = "") +
  coord_cartesian(ylim = c(-1,1))

# take mean from lm predictions ----

errors(preds_lm$dt_lm) <- preds_lm$se
mean_lm <- preds_lm %>%
  group_by(year) %>%
  summarise(mean_dt = log10(gm_mean(10^dt_lm)), .groups = "keep")
mean_lm$se <- errors(mean_lm$mean_dt) 
mean_lm$mean_dt <- drop_errors(mean_lm$mean_dt)

# plot!
ggplot(mean_lm, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_dt - 1.96*se,
                  ymax = mean_dt + 1.96*se)) +
  geom_line(aes(y = mean_dt), col = "white") +
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  labs(x = "")

#### hierarchical model (treating populations as covarying) ----

# build model
length(unique(ex_dt$year))/2
ex_dt$Binomial <- as.factor(ex_dt$Binomial)
hlm <- lm(size ~ 
              s(year, k = 13, bs = "tp") + 
              s(Binomial, k = 13, bs = "re"),
            data = ex_dt, method = "REML", family = "gaussian")

# extract covariance matrix
hlm_cov <- vcov(hlm, freq = TRUE) %>% cov2cor()
hlm_cov_pops <- hlm_cov[grep("Binomial", rownames(hlm_cov)),
                          grep("Binomial", colnames(hlm_cov))]
colnames(hlm_cov_pops) <- gsub("_", "\n", levels(ex_dt$Binomial))
rownames(hlm_cov_pops) <- gsub("_", "\n", levels(ex_dt$Binomial))
corrplot::corrplot(hlm_cov_pops, method = "color",
                   #addCoef.col = "white", # Add coefficient of correlation
                   #addCoefasPercent = TRUE,
                   tl.col="black", tl.srt = 50, tl.offset = 1)
# predict lm
preds_lm <- mgcv::predict.lm(hlm, se.fit = TRUE)

# add to data frame
ex_dt <- ex_dt %>%
  mutate(
    fit = preds_lm$fit,
    cilo = preds_lm$fit - 1.96*preds_lm$se.fit,
    cihi = preds_lm$fit + 1.96*preds_lm$se.fit
  )

# plot hlm trend
ggplot(ex_dt, aes(x = year)) +
  geom_ribbon(aes(ymin = cilo,
                  ymax = cihi)) +
  geom_line(aes(y = fit), col = "white") +
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  labs(x = "")