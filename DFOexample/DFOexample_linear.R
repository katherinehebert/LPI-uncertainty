# Script to demonstrate HMSC vs. a linear model on DFO dataset 
# for a Box demonstration of the benefits of hierarchical approaches

# clear workspace
rm(list=ls())

# set working directory
setwd("~/Documents/GitHub/LPI-sensitivity")
source('~/Documents/GitHub/LPI-sensitivity/scripts/scenario_functions.R')
# manually rewrite the hmsc predict function with script from a pull request
# on the repo, which adds quantiles to the predicted values
source('~/Documents/GitHub/LPI-sensitivity/DFOexample/hmsc_predict_PR27.R')

# load packages
library(tidyverse)
library(HMSC)
library(errors)
library(circlize)
library(DescTools)

# set ggplot theme
theme_set(theme_classic())

#### SET UP #### ---------------------------------------------------------------

# make vector of years we want to predict on
years <- data.frame("year" = c(1995:2016))

# read LPD dataset
lpd <- read_csv("data_raw/LPR2020data_public.csv")
ex_sub <- filter(lpd, Reference == "{DFO, 2016 #4093}", 
                 Location == "NAFO 4VW")

# clean up and convert to long format
ex <- ex_sub %>%
  pivot_longer(cols = c("1950":"2018"), 
               values_to = "size", 
               names_to = "year")
ex$year <- as.integer(ex$year)
ex$size <- as.numeric(ex$size)
#ex$size <- log10(ex$size)
ex <- drop_na(ex, size) # remove rows with NA
# subset to years with pretty complete data across many species
ex <- filter(ex, year %in% c(1995:2018))
ex <- filter(ex, ID != 22123)
# order ex according to year
ex <- ex[order(ex$year),]

# create clean wide version
ex_wide_og <- subset(ex, select = c(ID, year, size)) %>%
  pivot_wider(names_from = ID,
              values_from = size)
rownames(ex_wide_og) <- ex_wide_og$year  
ex_wide <- ex_wide_og[,-1]  
# remove columns with at least one NA
ex_wide <- ex_wide[ , colSums(is.na(ex_wide)) == 0]

#plot time series ----
(p_ts <- ggplot(ex) +
   geom_line(aes(x = year, y = log10(size), col = Binomial, group = ID), lwd = 0.3) +
   labs(y = "Abundance", col = "Population") +
   theme(legend.position = "none")
)

### GROWTH RATES ####

# function to calculate population growth rate (dt) using the chain method
# N = vector of population sizes (not transformed)
calc_dt_chain <- function(N){
  dt = c(0) # initial value
  for(i in 2:length(N)){
    dt[i] = log10(N[i]/N[i-1])
  }
  return(dt)
}

chain_dt <- list()
# calculate growth rates by chain method
for(i in 1:ncol(ex_wide)){
  chain_dt[[i]] <- calc_dt_chain(as.matrix(ex_wide)[,i])
}
ex_wide$year <- ex_wide_og$year


### AVERAGED LINEAR MODEL #### -------------------------------------------------

# run basic linear model of dt ~ year on each population
ex_lm <- list()
ex_lm_pred <- list()
for(i in 1:length(chain_dt)) {
  ex_lm[[i]] <- lm(chain_dt[[i]] ~ year, data = ex_wide)
  # predict model over years of interest
  ex_lm_pred[[i]] <- cbind(years,
                           predict(ex_lm[[i]], 
                                   years, 
                                   se.fit = TRUE)
  )
}
names(ex_lm_pred) <- colnames(ex_wide)[-ncol(ex_wide)]

# bind together
ex_lm_all <- bind_rows(ex_lm_pred, .id = "ID")
ex_lm_all$ID <- as.numeric(ex_lm_all$ID)

# plot growth rate predictions
ggplot(ex_lm_all) +
  geom_line(aes(x = year, y = fit, group = ID, col = ID)) +
  theme(legend.position = "none") +
  labs(y = "Predicted growth rates (log10)", x = "")


## TAKE GEOMETRIC MEAN GROWTH RATE ----

# first, match the population IDs with species names
ex_lm_all <- right_join(unique(subset(ex, select = c(ID, Binomial))), ex_lm_all)

# split into one group per year
ex_lm_all$year <- factor(ex_lm_all$year)
ex_lm_years <- split(ex_lm_all, f = ex_lm_all$year)

# take geometric mean with confidence intervals

# classic way
lm_means <- lapply(ex_lm_years, function(x) log10(Gmean(10^x$fit, conf.level = .95)))
lm_means <- bind_rows(lm_means, .id = "year")
# via bootstrapping
lm_boot <- lapply(ex_lm_years, function(x) log10(Gmean(10^x$fit, 
                                                       conf.level = .95, 
                                                       method = "boot")))
lm_boot <- bind_rows(lm_boot, .id = "year")

# error propagation source if it needs to be explained or done by hand:
# https://math.stackexchange.com/questions/123276/error-propagation-on-weighted-mean
gm_error <- function(gm, n, error_vec, value_vec){
  err <- abs((gm/n)*sqrt(sum((error_vec/value_vec)^2, na.rm = TRUE)))
}

# geometric mean per year
lm_df <- data.frame(
  "year" = years$year,
  "mean_dt" = NA,
  "error_dt" = NA
)
for(i in 1:nrow(lm_df)){
  lm_df$mean_dt[i] <- log10(Gmean(10^ex_lm_years[[i]]$fit))
  lm_df$error_dt[i] <- gm_error(gm = lm_df$mean_dt[i],
                                n = nrow(ex_lm_years[[i]]),
                                error_vec = ex_lm_years[[i]]$se.fit,
                                value_vec = ex_lm_years[[i]]$fit)
}

# plot the mean growth rate trend with error confidence interval 
ggplot(lm_df, aes(x = year)) +
  geom_ribbon(aes(ymin = mean_dt - 1.96*(error_dt), ymax = mean_dt + 1.96*(error_dt)), alpha = .3) +
  geom_line(aes(y = mean_dt), col = "white") + 
  labs(y = "Growth rate (log10)", x = "")

# this is the LPI version!
ggplot(lm_boot, aes(x = as.numeric(year))) +
  geom_ribbon(aes(ymin = lwr.ci, ymax = upr.ci)) +
  geom_line(aes(y = mean))

#### HMSC #### -----------------------------------------------------------------

#### Changed ex_l to just making a list called chain_dt. Adjust this part! ----

# make one data frame (which includes a column of growth rates)


# create year x species matrix  
Y <- bind_cols(chain_dt) %>% as.matrix()
rownames(Y) <- ex_wide$year
colnames(Y) <- colnames(ex_wide)[-ncol(ex_wide)]

# create XData dataframe of years
XData <- data.frame(x = ex_wide$year) 

# get number of populations & length of time series
npops <- ncol(Y)
tsl <- nrow(XData)


#### build HMSC model ####

# format data
XData <- as.matrix(XData)
YData <- as.matrix(Y)
randomEff <- factor(XData) %>% as.data.frame() 

# prepare data
formDat <- as.HMSCdata(X = XData, Y = YData, Random = randomEff)

# run model
m_hmsc <- hmsc(formDat,
                family = "gaussian",
                niter = 15000, #15000,
                nburn = 500, #5000,
                thin = 5,
                verbose = FALSE)
# save model
save(m_hmsc, file = "~/Documents/GitHub/LPI-sensitivity/DFOexample/hmsc.rds")

# extract the parameters and format them into an `mcmc` object.
mcmcMeansParamX <- as.mcmc(m_hmsc, parameters = "meansParamX")
# check for parameter convergence (visually)
traceplot(mcmcMeansParamX[,1])
traceplot(mcmcMeansParamX[,2])

# model explanatory power
(R2 <- Rsquared(m_hmsc, averageSp = FALSE))
(R2comm <- Rsquared(m_hmsc, averageSp = TRUE))

# average model
avg <- as.data.frame(apply(m_hmsc$results$estimation$paramX, 1:2, mean))
# get 95% confidence intervals
ci_lo <- as.data.frame(apply(m_hmsc$results$estimation$paramX, 1:2, quantile, probs = 0.025))
ci_hi <- as.data.frame(apply(m_hmsc$results$estimation$paramX, 1:2, quantile, probs = 0.975))

# # undo log transformation
# for(i in c(avg, ci_lo, ci_hi)){
#   i <- 10^i
# }

# plot all populations' models
ggplot(data = as.data.frame(XData), aes(x = V1)) +
  geom_abline(data = avg, 
              aes(slope = x, intercept = Intercept),
              lwd = 0.2) +
  geom_abline(data = ci_lo,
              aes(slope = x, intercept = Intercept),
              lty = 2, col = "blue", lwd = 0.1) +
  geom_abline(data = ci_hi,
              aes(slope = x, intercept = Intercept),
              lty = 2, col = "red", lwd = 0.1) +
  ylim(c(-1,1)) +
  theme_classic()

# get the LPI from the slopes of all populations at once -----------------------

# average model
avg_slope <- log10(Gmean(10^m_hmsc$results$estimation$meansParamX))
median_slope <- (median(m_hmsc$results$estimation$meansParamX))

# get 95% confidence intervals
ci_lo_slope <- quantile(m_hmsc$results$estimation$meansParamX, probs = 0.025)
ci_hi_slope <- quantile(m_hmsc$results$estimation$meansParamX, probs = 0.975)

# compute LPI with these averages
hmsc_df <- data.frame("mean_dt" = rep(avg_slope, nrow(XData)),
                     "cilo" = rep(ci_lo_slope, nrow(XData)),
                     "cihi" = rep(ci_hi_slope, nrow(XData)),
                     "year" = XData[,1])

#### species associations ####

# extract all estimated species-to-species associations matrix
assoMat <- corRandomEff(m_hmsc, cor = TRUE)

# Average
siteMean <- apply(assoMat[, , , 1], 1:2, mean)
# plot all of them!
#circlize::chordDiagram(siteMean, symmetric = TRUE)

# subset to associations that are >= 70%
siteMean_sub <- as.matrix(siteMean)
siteMean_sub[which(abs(siteMean_sub) < 0.7)] <- 0
# swap out the IDs for the binomial names
# colnames(siteMean_sub) <- unique(ex$Binomial[which(ex$ID %in% colnames(siteMean_sub))])
# rownames(siteMean_sub) <- unique(ex$Binomial[which(ex$ID %in% rownames(siteMean_sub))])

# plot styling
#pal = rainbow(ncol(siteMean_sub))
pal = viridis::viridis(ncol(siteMean_sub))
#pal = viridis::magma(ncol(siteMean_sub))
circlize::chordDiagram(siteMean_sub, symmetric = TRUE, 
                       annotationTrack = "grid", preAllocateTracks = 1, 
                       grid.col = pal)
circos.trackPlotRegion(track.index = 1, 
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         sector.name = get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), ylim[1] + .1, 
                                     sector.name, 
                                     facing = "clockwise", 
                                     niceFacing = TRUE, 
                                     adj = c(0, 0.5), 
                                     col = "black", 
                                     cex = .8)}, 
                       bg.border = NA)    
circos.clear()


# build matrix of colours for chordDiagram (red : neg, blue: pos)
siteDrawCol <- matrix(NA, nrow = nrow(siteMean), ncol = ncol(siteMean))
siteDrawCol[which(siteMean < 0, arr.ind=TRUE)]<-"red"
siteDrawCol[which(siteMean > 0, arr.ind=TRUE)]<-"blue"

# plot correlations (converted variance-covariances!!!)
corrplot::corrplot(siteMean, 
                   method = "color",
                   col = rev(colorRampPalette(c("blue","white","red"))(200)),
                   mar = c(0,0,1,0), tl.col = "grey40", tl.cex = 0.5)


#### predict and plot dt over time ####

# prepare and format data
formPredDat = as.HMSCdata(X = XData, Random = randomEff)

# predict values
pred <- predict(m_hmsc,
                newdata = formPredDat,
                type = "response")

# plot predicted values (one line per population)
# then plot mean dt trend across populations 
Colours <- rainbow(ncol(pred$mu))
par(mar = c(5.1, 4.1, 4.1, 2.1))
pred_mean = log10(apply(10^pred$mu, 1, gm_mean)) # take the mean dt per year
pred_cilo = log10(apply(10^pred$q025 , 1, gm_mean)) 
pred_cihi = log10(apply(10^pred$q975, 1, gm_mean)) 
plot(pred_mean ~ XData,
     col = "white", xlab = "years", ylab = "dt",
     ylim = c(min(pred$mu), max(pred$mu)))
for(i in 1:ncol(pred$mu)) {lines(pred$mu[,i] ~ XData,
                              col = Colours[i], lwd = 0.1)}
lines(pred_mean ~ XData, lwd = 1.2)
lines(pred_cilo ~ XData, lty = 2)
lines(pred_cihi ~ XData, lty = 2)
for(i in 1:ncol(pred$mu)) {lines(chain_dt[[i]] ~ XData,
                                 col = Colours[i], lwd = 0.2)}

## comparison plot of linear vs. hmsc ####


# compare predictions here

pred_lm <- data.frame(
  "year" = lm_boot$year,
  "mean_dt" = lm_boot$mean,
  "cilo" = lm_boot$lwr.ci,
  "cihi" = lm_boot$upr.ci,
  "model" = "lm"
)
pred_hmsc <- data.frame(
  "year" = XData[,1],
  "mean_dt" = pred_mean,
  "cilo" = pred_cilo,
  "cihi" = pred_cihi,
  "model" = "HMSC"
)
# combine
pred_results <- rbind(pred_lm, pred_hmsc)
# plot
ggplot(pred_results, aes(x = year, group = model)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = model), 
              lwd = 0, alpha = .3) +
  geom_line(aes(y = mean_dt, col = model)) +
  ylim(-1,1)
# seems like there's a log10 issue here. check. ----

# compare the slopes themselves here (with their associated errors)
# simple slope and intercept line comparisons

# take mean slope of linear models
lm_res <- lapply(ex_lm, coef) %>% bind_rows()
lm_res <- Gmean(10^lm_res$year, conf.level = .95) %>% log10() %>% t() %>% as.data.frame()
hmsc_res <- hmsc_df[1,1:3]
colnames(hmsc_res) <- colnames(lm_res)

# this is all in log10
hmsc_res$model <- "HMSC"
lm_res$model <- "lm"
# bind together
slope_res <- rbind(lm_res, hmsc_res)

# boxplot to compare them
ggplot(slope_res) +
  geom_point(aes(x = model, y = mean)) +
  geom_errorbar(aes(x = model, ymin = lwr.ci, ymax = upr.ci), 
                width = .1, lwd = .2) +
  ylim(c(-1,1))

# residual error
m_hmsc$results$estimation$varX
# get residaul variance from lms

resid(ex_lm[[1]])
resid(m_hmsc)
