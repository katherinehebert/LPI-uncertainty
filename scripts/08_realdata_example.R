# real data example to demonstrate the mathematical analysis of error 
# propagation in the LPI

# set-up ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

theme_set(ggpubr::theme_pubr())

# read LPD dataset
lpd <- read_csv("data_raw/LPR2020data_public.csv")

# subset to example system
ex_sub <- lpd %>% 
  #filter(Location == "Kluane Lake Research Station, Yukon, Canada") 
  filter(ID %in% c(4363, 4362))

# clean up and convert to long format
ex <- ex_sub %>%
  pivot_longer(cols = c("1950":"2018"), 
               values_to = "N", 
               names_to = "year")
ex$year <- as.integer(ex$year)
ex$N <- as.numeric(ex$N)
ex <- drop_na(ex, N) # remove rows with NA
ex_l <- ex

# create clean wide version
ex <- ex_sub[,30:ncol(ex_sub)] %>% 
  apply(1:2, as.numeric) %>%
  t()
colnames(ex) = ex_sub$Binomial
ex <- na.omit(ex)

# EQ. 1 get step-wise growth rates ----

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
dt <- as.data.frame(ex)
for(i in 1:2){
  dt[,i] <- calc_dt_chain(ex[,i])
}

# set up to plot growth rates 
dt$year <- as.integer(rownames(ex))
ex_dt <- as.data.frame(dt) %>%
  pivot_longer(cols = 1:2, names_to = "Binomial", values_to = "size")

## Uncertainty propagation ----

# EQ 2. Expected growth rate ----
# expectation of the growth rate trend
eq2 <- function(N, sigma_measure = 0){
  d <- c()
  for(t in 2:length(N)){
    d[t] <- log10(N[t]/N[t-1]) + (sigma_measure^2)/(2*(N[t-1]^2 - N[t]^2))
  }
  return(d)
}

sigma1 = var(log10(ex[,1]))
sigma2 = var(log10(ex[,2]))
              
dt <- cbind(eq2(ex[,1], sigma_measure = sigma1),
            eq2(ex[,2], sigma_measure = sigma2))

# equation 2
# just the growth rate part
eq2_donly <- function(N){
  eq2dt <- c()
  for(t in 2:length(N)){
    eq2dt[t] <- log10(N[t]/N[t-1])
  }
  return(eq2dt)
}
donly <- cbind(eq2_donly(ex[,1]), eq2_donly(ex[,2]))

# equation 2
# just the uncertainty part
eq2_varonly <- function(N, sigma_measure = 0){
  eq2var <- c()
  for(t in 2:length(N)){
    eq2var[t] <- (sigma_measure^2)/(2*((N[t-1])^2 - (N[t])^2))
  }
  return(eq2var)
}
varonly <- cbind(eq2_varonly(ex[,1], sigma_measure = sigma1),
                 eq2_varonly(ex[,2], sigma_measure = sigma2))


## EQ 3 - flag ----
## we do not have process uncertainty. just assuming everything is measurement uncertainty

eq3 <- function(N, sigma_measure = 0, sigma_process = 0){
  eq3var <- c()
  for(t in 2:length(N)){
    eq3var[t] <- sigma_process^2 + (sigma_measure^2)*((1/((N[t])^2) - (1/(N[t-1])^2)))
  }
  return(eq3var)
}
var_dt <- cbind(eq3(ex[,1], sigma_measure = sigma1, sigma_process = 0),
                eq3(ex[,2], sigma_measure = sigma2, sigma_process = 0))

## EQ. 4 -----
## average growth rate trend
dt_bar <- apply(dt, 1, mean)
dt_bar[1] = 0

## EQ. 5 ----
## uncertainty in the average growth rate

dt_cov <- cov(dt, use = "pairwise.complete.obs")
dt_cov <- dt_cov[which(lower.tri(dt_cov))]
var_dtbar = (1/nrow(ex))*(apply(var_dt, 1, sum) + 2*sum(dt_cov))
var_dtbar[1] = 0

## EQ. 6 ----
## computation of the LPI

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
var_lpi[1] = 0
plot(var_lpi, type = "l")


# results to save

lpi_res <- data.frame(
  "time" = as.integer(rownames(ex)),
  "dtbar" = dt_bar,
  "dtbar_variance" = var_dtbar,
  "lpi_nocorrection" = lpi_nocorrection, # eq 6
  "lpi_correction" = lpi_correction, # eq 7
  "lpi_variance" = var_lpi, # eq 8,
  "lpi_bias" = lpi_nocorrection - lpi_correction
)
saveRDS(lpi_res, "outputs/realdata_example_uncertaintypropagation.RDS")


## plots ----

ex_dt$Binomial <- gsub("_", " ", ex_dt$Binomial)
ex_l$Binomial <- gsub("_", " ", ex_l$Binomial)

# plot time series
(p_ts <- ggplot(ex_l) +
   geom_line(aes(x = year, y = N, col = Binomial), lwd = .8) +
   labs(y = "Abundance", x = "", col = "Population", title = "Raw data") +
   theme(legend.position = "none")
)

(p_dt <- ggplot(ex_dt) +
  geom_line(aes(x = year, y = size, col = Binomial), lty = 2, lwd = .5) +
  geom_line(data = lpi_res, aes(x = time, y = dtbar)) +
  geom_hline(yintercept = 0, lwd = .1) +
  labs(y = "Average growth rate", x = "", col = "Population", title = "Equation 2 & 4") +
    theme(legend.position = "none"))

(p_dtbar <- ggplot(ex_dt) +
    geom_line(data = lpi_res, aes(x = time, y = dtbar_variance)) +
    labs(y = "Variance in the \naverage growth rate", x = "", title = "Equation 5"))

(LPI = (ggplot(lpi_res) +
              geom_line(aes(x = time, y = lpi_nocorrection), lty = 3, lwd = .4) +
              geom_line(aes(x = time, y = lpi_correction)) +
              labs(y = "LPI", x = "", title = "Equation 6 & 7") 
))
(bias = (ggplot(lpi_res) +
   geom_line(aes(x = time, y = lpi_bias)) +
     geom_hline(yintercept = 0, lwd = .1) +
     labs(y = "Bias in the LPI", x = "", title = "Equation 7") 
))
(variance = (ggplot(lpi_res) +
    geom_line(aes(x = time, y = lpi_variance)) +
    labs(y = "Variance of the LPI", x = "", title = "Equation 8") 
))

library(patchwork)
(p_ts + p_dt + p_dtbar)/ (LPI + bias + variance) + plot_annotation(tag_levels = "a")
ggsave("figures/realdata_example_uncertaintypropagation.png", width = 10, height = 6.5)
