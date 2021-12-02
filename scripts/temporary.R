
library(ggridges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)

# compare distribution of abundances

# N with observation error
Nerror <- readRDS("~/Documents/GitHub/LPI-sensitivity/simulations/scenario2A_l.RDS")
Nerror <- readRDS("~/Documents/GitHub/LPI-sensitivity/simulations/scenario2B_l.RDS")
Nerror <- readRDS("~/Documents/GitHub/LPI-sensitivity/simulations/scenario2C_l.RDS")

# N without observation error
Ntrue <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/scenario2A_true.RDS")
Ntrue <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/scenario2B_true.RDS")
Ntrue <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/scenario2C_true.RDS")

colnames(Ntrue)[5] <- "Ntrue"
colnames(Nerror)[5] <- "Nerror"

N <- inner_join(Ntrue, Nerror) %>%
  pivot_longer(cols = c("Ntrue", "Nerror"), 
               values_to = "N", names_to = "N_method")

dlnorm(quantile(Ntrue$Ntrue))
dlnorm(quantile(Nerror$Nerror))

sd(Nerror$Nerror)/mean(Nerror$Nerror)
sd(Ntrue$Ntrue)/mean(Ntrue$Ntrue)


ggplot(N, aes(x = N, y = N_method)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975), alpha = 0.7,
                      jittered_points = TRUE, position = "raincloud")

## compare distribution of growth rates

# from rlpi calculation on ERROR N
dt_error <- read_csv("outputs/rlpi/scenario1C_pops_lambda.csv", col_types = cols(X1 = col_skip())) %>%
  pivot_longer(cols = 3:13, values_to = "dt_error", names_to = "time")

# from rlpi calculation on TRUE N
dt_true <- read_csv("outputs/rlpi/scenario1C_true_pops_lambda.csv",  col_types = cols(X1 = col_skip())) %>% 
  pivot_longer(cols = 3:13, values_to = "dt_true", names_to = "time")

# from chain calculation on raw abundances
dt_chain <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/scenario1C_l_dtchain.RDS")
dt_chain$time <- (dt_chain$time + 1969) %>% as.character()

ggplot(filter(Nerror, time != 1), aes(x = Nerror, y = factor(time))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.5), alpha = 0.7)

ggplot(filter(dt_chain, time != 1), aes(x = dt_chain, y = factor(time))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.5), alpha = 0.7) +
  coord_cartesian(xlim = c(-0.2, 0.2))
ggplot(filter(dt_error, time != "1970"), aes(x = dt_error, y = factor(time))) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.5), alpha = 0.7) +
  coord_cartesian(xlim = c(-0.2, 0.2))

# compare growth rates between methods
dt <- inner_join(dt_true, dt_error) %>% 
  inner_join(dt_chain) %>%
  pivot_longer(cols = c(dt_true, dt_error, dt_chain), 
               values_to = "dt", names_to = "method")

ggplot(filter(dt, time != "1970"), aes(x = 10^dt, y = method)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.5), alpha = 0.7,
                      jittered_points = TRUE, position = "raincloud", scale = 0.5) +
  geom_vline(data= filter(dt_true, time != "1970"), aes(xintercept = EnvStats::geoMean(10^dt_true)), col = "red") 

ggplot(filter(dt, time != "1970")) +
  geom_boxplot(aes(x = method, y = dt, fill = method)) +
  theme_linedraw() 

compare_dists <- function(simID){
  
  # from rlpi calculation on ERROR N
  dt_error <- read_csv(paste0("outputs/rlpi/scenario", simID, "_pops_lambda.csv"), col_types = cols(X1 = col_skip())) %>%
    pivot_longer(cols = 3:13, values_to = "dt_error", names_to = "time")
  
  # from rlpi calculation on TRUE N
  dt_true <- read_csv(paste0("outputs/rlpi/scenario", simID, "_true_pops_lambda.csv"),  col_types = cols(X1 = col_skip())) %>% 
    pivot_longer(cols = 3:13, values_to = "dt_true", names_to = "time")
  
  # from chain calculation on raw abundances
  dt_chain <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/scenario", simID, "_l_dtchain.RDS"))
  dt_chain$time <- (dt_chain$time + 1969) %>% as.character()
  
  # compare growth rates between methods
  dt <- inner_join(dt_true, dt_error) %>% 
    inner_join(dt_chain) %>%
    pivot_longer(cols = c(dt_true, dt_error, dt_chain), 
                 values_to = "dt", names_to = "method")
  
  dt1 = 10^(dt_true$dt_true[-which(dt_true$time == "1970")])
  dt2 = 10^(dt_error$dt_error[-which(dt_error$time == "1970")]) 
  dt3 = 10^(dt_chain$dt_chain[-which(dt_chain$time == "1970")])
  
  # dlnorm(quantile(dt1))
  # dlnorm(quantile(dt2))
  # dlnorm(quantile(dt3))
  # 
  # sd(dt1)
  # sd(dt2)
  # sd(dt3)
  
  samp1 = dt3
  samp2 = dt2
  
  n = length(samp1)
  m = length(samp2)
  
  # Fn vs. Gm
  plot(ecdf(samp1), ylab = "Probability", xlim = c(0.8,1.3), 
       xlab = "Growth rate (dt)", main = paste("scenario", simID))
  lines(ecdf(samp2), main = "", col = 4)
  
  # Add Dnm1
  samp1_sorted <- sort(samp1)
  samp2_sorted <- sort(samp2)
  Dnm_1 <- abs((1:n) / n - ecdf(samp2)(samp1_sorted))
  i1 <- which.max(Dnm_1)
  lines(rep(samp2_sorted[i1], 2), 
        c(i1 / m, ecdf(samp1_sorted)(samp2_sorted[i1])), 
        col = 3, lwd = 2, type = "o", pch = 16, cex = 0.75)
  rug(samp1, col = 1)
  
  # Add Dnm2
  Dnm_2 <- abs(ecdf(samp1)(samp2_sorted) - (1:m) / m)
  i2 <- which.max(Dnm_2)
  lines(rep(samp1_sorted[i2], 2), 
        c(i2 / n, ecdf(samp2_sorted)(samp1_sorted[i2])), 
        col = 2, lwd = 2, type = "o", pch = 16, cex = 0.75)
  rug(samp2, col = 4)
  
  # are the distributions different?
  ks.test(x = samp1, y = samp2) 
  
}

triplot <- function(scenario_number){
  par(mfrow = c(1,4))
  compare_dists(paste0(scenario_number, "A"))
  compare_dists(paste0(scenario_number, "B"))
  compare_dists(paste0(scenario_number,"C"))
  # add legend plot
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("topleft", lwd = 2, col = c(1, 4, 3, 2), 
         legend =
           latex2exp::TeX(c("$dt_{log-ratio}$", "$dt_{rlpi}$", "$D_{log-ratio,rlpi,1}$", "$D_{log-ratio,rlpi,2}$")))
}

# plot distribution comparisons in one pdf file
pdf(file = here("outputs/dt_KSComparisons.pdf"))
  sapply(as.character(1:7), FUN = triplot)
dev.off()

triplot("1")
triplot("2")
triplot("3")

hist(Nerror$Nerror)
sd(Nerror$Nerror)/mean(Nerror$Nerror)
x <- filter(dt_error, dt_error < 1) 
hist(x$dt_error)
10^sd(x$dt_error)

mean(Nerror$Nerror)
mean()
# testing geometric mean
# there's a bigger gap when there's more variance..... but like not huge
# although could try with 10^ then logging, but hm.

library(EnvStats)
n = 100
sds <- seq(from = 0.05, to = 0.25, by = 0.05)
ex <- matrix(NA, nrow = n, ncol = length(sds))

par(mfrow = c(3, 5))
for(i in 1:length(sds)){
  ex[,i] <- rnorm(n, mean = 0.7, sds[i])
  hist(ex[,i], xlim = c(0, 2), main = paste("SD = ", sds[i]))
  abline(v = geoMean(ex[,i]), col = "red")
  abline(v = mean(ex[,i]), col = "blue")
}
for(i in 1:length(sds)){
  ex[,i] <- rnorm(n, mean = 1, sds[i])
  hist(ex[,i], xlim = c(0, 2), main = paste("SD = ", sds[i]))
  abline(v = geoMean(ex[,i]), col = "red")
  abline(v = mean(ex[,i]), col = "blue")
}
for(i in 1:length(sds)){
  ex[,i] <- rnorm(n, mean = 1.25, sds[i])
  hist(ex[,i], xlim = c(0, 2), main = paste("SD = ", sds[i]))
  abline(v = geoMean(ex[,i]), col = "red")
  abline(v = mean(ex[,i]), col = "blue")
}



mean(dt_chain$dt_chain)
sd(dt_chain$dt_chain)

mean(dt_error$dt_error)
sd(dt_error$dt_error)


# does the distribution of rlpi dts match the distribution of chain method dts?
# (where chain method is the raw error)
dt_precisiontest <- dt %>% filter(!method == "dt_true")
t.test(dt ~ method, data = dt_precisiontest, var.equal = F)

# get sd of the rlpi
rlpi <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/rlpi/scenario1C_rlpi.RDS")

(((rlpi$CI_high-rlpi$CI_low)/3.92)*sqrt(2000)) %>% mean()
(rlpi$CI_high - rlpi$LPI_final)/1.96

# try when wifi is possible
library(brms)
library(broom.mixed)

CHAINS <- 4
ITER <- 2000
WARMUP <- 1000
BAYES_SEED <- 1234
options(mc.cores = parallel::detectCores())  # Use all cores

# model the relationship between dt and method, and the variance and method

# think of priors
par(mfrow=c(1,1))
curve(expr = dcauchy(x, location = 0.05, scale = .1), from = -5, to = 5)


# calculate the medians of the posterior distributions and create confidence intervals

brms_uneq <- brm(
  bf(dt ~ method, sigma ~ method), 
  data = mutate(dt_precisiontest, method = fct_rev(method)),
  prior = c(set_prior("normal(0, 2)", class = "Intercept"),
            set_prior("normal(0, 2)", class = "b"),
            set_prior("cauchy(0.05, .1)", class = "b", dpar = "sigma")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  file = "cache/brms_uneq"
)

brms_uneq_tidy <- 
  tidyMCMC(brms_uneq, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")
brms_uneq_tidy
# intercept = mean dt score
# 95% certain that the difference in dt between methods is between -0.012 and 0.011,
# with a median of -0.00005. In other words, the difference is tiny


