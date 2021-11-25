
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
  dt_compare <- inner_join(bootstrap, chain, by = "year") %>%
    pivot_longer(cols = c(bootstrap, chain), names_to = "type", values_to = "lpi")
  
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
  df$true_diff_rlpiCI_low <- NA
  df$true_diff_rlpiCI_high <- NA
  for(i in 1:nrow(df)){
    if(df$true_within_rlpiCI[i] == "Failure"){
      if(df$LPI_final_true[i] < df$q025_rlpi[i]){
        df$true_diff_rlpiCI_low[i] <- df$q025_rlpi[i] - df$LPI_final_true[i]
      } else if(df$LPI_final_true[i] > df$q975_rlpi[i]) {
        df$true_diff_rlpiCI_high[i] <- df$LPI_final_true[i] - df$q975_rlpi[i] 
      }
    } else next
  }
  
  # calculate interval widths & difference in interval width
  df$rlpiCI_width <- df$q975_rlpi - df$q025_rlpi
  df$chainCI_width <- df$q975_chain - df$q025_chain
  df$CI_diff <- df$chainCI_width - df$rlpiCI_width
  saveRDS(df, paste0("outputs/scenario", simID, "_precision.RDS"))
  
}
compare_precision("1A")
compare_precision("1B")
compare_precision("1C")
compare_precision("2A")
compare_precision("2B")
compare_precision("4A")

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
lapply(sim_ids, compare_precision)
