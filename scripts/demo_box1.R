# Script to demonstrate the need for a multispecies index of biodiversity change
# using a simple ecological example

# set-up ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(hrbrthemes)
library(patchwork)
library(synchrony)
library(errors)
library(mgcv)
library(gratia)

# load required functions
source('scripts/scenario_functions.R')

# function to calculate population growth rate (dt) using the chain method
# N = vector of population sizes (not transformed)
calc_dt_chain <- function(N){
  dt = c(0) # initial value
  for(i in 2:length(N)){
    dt[i] = log10(N[i]/N[i-1])
  }
  return(dt)
}
# set ggplot theme
theme_set(ggpubr::theme_pubr() +
            theme(legend.position = c(0.8, 0.92),
                  legend.direction = "vertical",
                  axis.title = element_text(size = 18),
                  axis.text = element_text(size = 16),
                  legend.text = element_text(size = 11))) 




# 1. explore data to find a good example ----

# read LPD dataset
lpd <- read_csv("data_raw/LPR2020data_public.csv")

# subset to example system
ex_sub <- lpd %>% 
  filter(Location == "Kluane Lake Research Station, Yukon, Canada") 

# clean up and convert to long format
ex <- ex_sub %>%
  pivot_longer(cols = "1950":ncol(.), 
               values_to = "size", 
               names_to = "year")
ex$year <- as.integer(ex$year)
ex$size <- as.numeric(ex$size)
ex <- drop_na(ex, size) # remove rows with NA
  
# create clean wide version
ex_wide <- ex_sub[,30:ncol(ex_sub)] %>% 
  apply(1:2, as.numeric) %>%
  t()
colnames(ex_wide) = ex_sub$Binomial

# plot time series ----
(p_ts <- ggplot(ex) +
  geom_line(aes(x = year, y = size, col = Binomial), lwd = 1.1) +
  scale_y_log10() +  
  coord_cartesian(xlim = c(1990, max(ex$year))) + # crop 
   scale_color_brewer(palette = "Dark2") +
  theme_ipsum_rc() +
  labs(y = "Population size (log10)"
       #title = "Interacting populations",
       #subtitle = "Coyote, Snowshoe hare, and Canadian lynx"
       )
)

# get step-wise growth rates ----

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
  scale_color_brewer(palette = "Dark2") +
  theme_ipsum_rc() +
  labs(y = "Growth rate (log10)"
       #title = "Interacting populations",
       #subtitle = "Coyote, Snowshoe hare, and Canadian lynx"
  )

(p_ts / p_dt)


# plot dt(t) ~ d(t+1) - stability plot ----------------------------------------

dt_stability <- ex_dt %>%
  group_by(Binomial) %>%
  summarise(
    year = year,
    dt_0 = size,
    dt_1 = size,
    )

p_all <- ggplot(data = dt_stability) +
  geom_point(aes(x = dt_0, y = lag(dt_1), col = Binomial)) +
  geom_abline(aes(slope = 1, intercept = 0), lwd = .2)  +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")

p_each <- p_all + facet_wrap(~Binomial)

(p_all / p_each)

# calculate linear covariance between populations ----

cov_ex_dt <- dt[,2:4] %>% cov() %>% cov2cor()
corrplot::corrplot(cov_ex_dt, type = "upper", method = "color",
                   addCoef.col = "white", # Add coefficient of correlation
                   tl.col="black", tl.srt = 45)


# # calculate synchrony ----
# 
# trio_sync <- community.sync(dt[,2:4], nrands = 100)
# plot(trio_sync)
# 

#### mean growth rate per time step (chain) ----
# (not treating populations as covarying)

mean_dt <- data.frame(
  year = dt[,1],
  dt = log10(apply(10^dt[,2:4], 1, EnvStats::geoMean)),
  sd = log10(apply(10^dt[,2:4], 1, EnvStats::geoSD))
)

# plot mean growth rate 
(p_mean <- ggplot(mean_dt) +
  geom_ribbon(aes(x = year, 
                  ymin = dt - 1.96*sd,
                  ymax = dt + 1.96*sd)) +
  geom_line(aes(x = year, y = dt), color = "white", lwd = 1.1) +
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  #scale_y_log10() +  
  theme_ipsum_rc() +
  labs(y = "Growth rate (log10)")
)

# simple linear model per time series

df <- as.data.frame(dt)

gam_canis <- gam(Canis_latrans ~ s(year), 
                 data = df,
                 method = "REML")
draw(gam_canis, residuals = TRUE)
mgcv::gam.check(gam_canis)

gam_lepus <- gam(Lepus_americanus ~ s(year), 
                 data = df, 
                 method = "REML")
draw(gam_lepus, residuals = TRUE)
mgcv::gam.check(gam_lepus)

gam_lynx <- gam(Lynx_canadensis ~ s(year),
                 data = df, 
                 method = "REML")
draw(gam_lynx, residuals = TRUE)
mgcv::gam.check(gam_lynx)

# predict each gam individually
pred_canis <- predict(gam_canis, se = TRUE)
pred_lepus <- predict(gam_lepus, se = TRUE)
pred_lynx <- predict(gam_lynx, se = TRUE)
# organise into a data frame
preds_gam <- data.frame(
  population = rep(colnames(dt)[2:4], each = nrow(dt)),
  year = rep(dt[,1], 3),
  dt_gam = c(pred_canis$fit, pred_lepus$fit, pred_lynx$fit),
  se = c(pred_canis$se.fit, pred_lepus$se.fit, pred_lynx$se.fit)
)

# rename columns before merging with prediction data frame
colnames(ex_dt)[2] <- "population"
colnames(ex_dt)[3] <- "true_mean"
preds_gam <- full_join(preds_gam, ex_dt)
# plot gam predictions vs. "true" growth rates
(p_gam_each <- ggplot(data = preds_gam, aes(x = year)) +
    geom_line(aes(y = true_mean, col = population)) +
    geom_ribbon(aes(x = year,
                    ymin = dt_gam - 1.96*se,
                    ymax = dt_gam + 1.96*se,
                    fill = population),
                alpha = .5) +
    geom_line(aes(y = dt_gam), col = "white") +
    #geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
    facet_wrap(~population) +
    labs(x = "", y = "Growth rate (log10)") +
    #coord_cartesian(ylim = c(-.6,.6)) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_x_continuous(breaks = c(1990, 2000, 2010)) +
    guides(fill = FALSE) +
    theme_ipsum_rc(plot_margin = margin(5, 5, 5, 5)) +
    theme(legend.position = "none")
)


# take mean from gam predictions ----

errors(preds_gam$dt_gam) <- preds_gam$se
mean_gam <- preds_gam %>%
  group_by(year) %>%
  summarise(mean_dt = log10(gm_mean(10^dt_gam)), .groups = "keep")
mean_gam$se <- errors(mean_gam$mean_dt) 
mean_gam$mean_dt <- drop_errors(mean_gam$mean_dt)

# plot!
(p_gam_mean <- ggplot() +
    geom_line(data = mean_dt, 
              aes(x = year, y = dt)) +
    geom_ribbon(data = mean_gam, 
                aes(x = year,
                    ymin = mean_dt - 1.96*se,
                    ymax = mean_dt + 1.96*se),
                alpha = .5) +
    geom_line(data = mean_gam, 
              aes(x = year, y = mean_dt), 
              col = "white") +
    geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
    labs(x = "", y = "Growth rate (log10)") +
    coord_cartesian(ylim = c(-1.1, 1.1)) +
    theme_ipsum_rc(plot_margin = margin(5, 5, 5, 5))
)



#### hierarchical model (treating populations as covarying) ----

# build model
df <- dt_stability[,1:3]
df$Binomial <- as.factor(df$Binomial)

# model GI
hgam <- gam(dt_0 ~ s(year, bs = "tp") +  # global part
              s(year, by = Binomial, m = 1, bs = "tp") + # individual smoothers per population
              s(Binomial, bs = "re"),
            data = df, 
            #sp = c(1, 1, 1, 1, 100),
            method = "REML")
gratia::appraise(hgam)
plot(hgam, residuals = TRUE)
draw(hgam, select = 1:5, scales = "fixed")

# extract global smoother trend
hgam_globaltrend <- evaluate_smooth(hgam, "year", n = 20)
# plot the smoother by year
ggplot(data = hgam_globaltrend[1:20,], aes(x = year)) +
  geom_ribbon(aes(ymin = est -1.96*se,
                  ymax = est + 1.96*se)) +
  geom_line(aes(y = est))

## covariance insights ----

# extract covariance matrix
hgam_cov <- vcov(hgam, freq = TRUE) %>% cov2cor()
hgam_cov_pops <- hgam_cov[grep("Binomial", rownames(hgam_cov)),
                          grep("Binomial", colnames(hgam_cov))]
# colnames(hgam_cov_pops) <- gsub("Binomial", "", colnames(hgam_cov_pops))
# rownames(hgam_cov_pops) <- gsub("Binomial", ":", rownames(hgam_cov_pops))

# slopes by year (GLOBAL TREND)
corrplot::corrplot(hgam_cov[1:9,1:9], 
                   method = "color", 
                   diag = FALSE, 
                   type = "lower",
                   #addCoef.col = "white", # Add coefficient of correlation
                   #addCoefasPercent = TRUE,
                   tl.cex = 0.5, tl.col="black", 
                   #tl.srt = 45, 
                   tl.offset = 1)

# intercepts
corrplot::corrplot(hgam_cov[35:37,35:37], 
                   method = "color", 
                   diag = FALSE, 
                   type = "lower",
                   #addCoef.col = "white", # Add coefficient of correlation
                   #addCoefasPercent = TRUE,
                   tl.cex = 0.5, tl.col="black", 
                   #tl.srt = 45, 
                   tl.offset = 1)

# year smoothers per population
corrplot::corrplot(hgam_cov[10:34,10:34], 
                   method = "color", 
                   diag = FALSE, 
                   type = "lower",
                   #addCoef.col = "white", # Add coefficient of correlation
                   #addCoefasPercent = TRUE,
                   tl.cex = 0.5, tl.col="black", 
                   #tl.srt = 45, 
                   tl.offset = 1)

# predict gam
preds_hgam <- mgcv::predict.gam(hgam, se.fit = TRUE)
# extract global trend only
preds_gam_G <- data.frame(
  year = unique(df$year),
  fit = preds_hgam$fit[1:26],
  cilo = preds_hgam$fit[1:26] - 1.96*preds_hgam$se.fit[1:26],
  cihi = preds_hgam$fit[1:26] + 1.96*preds_hgam$se.fit[1:26]
)

# plot hgam trend
(p_hgam <- ggplot() +
  geom_line(data = mean_dt, 
            aes(x = year, y = dt)) +
  geom_ribbon(data = preds_gam_G, 
              aes(x = year,
                      ymin = cilo,
                  ymax = cihi),
              alpha = .5) +
  geom_line(data = preds_gam_G, 
            aes(x = year, 
                y = fit), 
            col = "white") +
  geom_hline(yintercept = 0, lty = 2, lwd = .2, col = "white") +
  labs(x = "", y = "Growth rate (log10)") +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
    scale_color_brewer(palette = "Dark2") +
  theme_ipsum_rc(plot_margin = margin(5, 5, 5, 5))
)


# independent species vs. multispecies approach --------------------------------

# (p_dt_lepus <-  ggplot(filter(ex_dt, Binomial == "Lepus_americanus")) +
#   geom_line(aes(x = year, y = size, col = Binomial), lwd = 1.1) +
#   #scale_y_log10() +  
#   coord_cartesian(xlim = c(1990, max(ex$year)),
#                   ylim = c(-1.5, 1.5)) + # crop 
#   ggpubr::theme_pubr() +
#   labs(y = "Growth rate (log10)", x = "",
#        col = "") +
#    theme(legend.position = c(0.8, 0.95),
#          legend.direction = "vertical",
#          axis.title = element_text(size = 18),
#          axis.text = element_text(size = 16),
#          legend.text = element_text(size = 11))
# )
# ggsave("figures/box1_ts_lepus.png", width = 5, height = 5)


(p_dt_web <-  ggplot(ex_dt) +
    geom_line(aes(x = year, y = true_mean, col = population), 
              lwd = 1.1) +
    coord_cartesian(xlim = c(1990, max(ex$year)),
                    ylim = c(-1.5, 1.5)) + # crop 
    ggpubr::theme_pubr() +
    labs(y = "Growth rate (log10)", x = "",
         col = "") +
    theme_ipsum_rc(plot_margin = margin(5, 5, 5, 5))+
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 11))
)
ggsave("figures/box1_ts_web.png", width = 5.7, height = 3.12)

p_dt_each <- p_dt_web + 
  facet_wrap(~population) + 
  scale_x_continuous(breaks = c(1990, 2000, 2010)) +
  guides(col = FALSE) + labs(title = "") +
  theme_ipsum_rc(plot_margin = margin(5, 5, 5, 5))
p_dt_mean <- p_mean + theme(plot.title = element_blank())


p_dt_each / p_gam_each / p_gam_mean + plot_annotation(tag_levels = "a")
#ggsave()


# comparison of mean gam and hgam
p_gam_mean + p_hgam + plot_annotation(tag_levels = "a")
