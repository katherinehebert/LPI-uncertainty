## checks: where does the error come from in the LPI calculation?

# load packages and prep environment
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(patchwork)
theme_set(ggpubr::theme_pubr())

# # import expected uncertainties of each scenario's growth rates
# uncertainties <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_uncertainty.RDS")), readRDS)
# names(uncertainties) <- gsub("_uncertainty.RDS", "", list.files(path = "outputs/", pattern = "_uncertainty.RDS"))
# uncertainty <- bind_rows(uncertainties, .id = "scenario") %>% subset(select = "scenario")

# import results
results <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_results.RDS")[-1]), readRDS)
names(results) <- gsub("_results.RDS", "", list.files(path = "outputs/", pattern = "_results.RDS")[-1])
results <- bind_rows(results, .id = "scenario")

# read results file
df <- readRDS("outputs/all_results.RDS")

# format interaction to show up on the plots in the right order
df$interaction <- gsub("-0.2", "Strong Synchrony", df$interaction) 
df$interaction <- gsub("-0.1", "Weak Synchrony", df$interaction) 
df$interaction <- gsub("0.1", "Weak Asynchrony", df$interaction) 
df$interaction <- gsub("0.2", "Strong Asynchrony", df$interaction) 
df$interaction <- gsub("0", "No Synchrony", df$interaction) 
df$interaction <- factor(df$interaction, 
                         levels = c("Strong Asynchrony", "Weak Asynchrony", "No Synchrony", "Weak Synchrony", "Strong Synchrony"))
# same with process error
df$Process_error <- factor(df$Process_error, levels = c("0", "0.1", "0.2"))

# combine with results including standard error for each step
results$scenario <- uncertainty$scenario
df_err <- inner_join(results, df)

## comparing GAM predictions to simulated population sizes

facet_names <- c(
  `0` = "Process ε = 0",
  `0.1` = "Process ε = 0.1",
  `0.2` = "Process ε = 0.2"
)

# plot difference between prediction and simulation
ggplot(filter(df_err, Lag == "0")) +
  geom_jitter(aes(y = (10^N_pred)-(10^N), 
                  x = interaction, 
                  col = direction), size = 1, alpha = .2, position = position_jitterdodge()) +
  facet_wrap(~ Process_error, dir = "v", labeller = as_labeller(facet_names)) +
  labs(x = "", y = expression(N[sim]~Delta~N[GAM]), col = "Trend") +
  theme(legend.position = "top") +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,3,5)]) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  ggpubr::theme_pubr() +
  geom_hline(yintercept = 0, lwd = .3, lty = 2)
ggsave("figures/figsupp_GAMpredictions.png", width = 6, height = 6)

# plot N from GAM prediction vs. simulated N
ggplot(filter(df_err, Lag == "0" & Process_error == "0")) +
  geom_point(aes(y = 10^N_pred, 
                  x = 10^N, 
                  col = direction), size = 1, alpha = .7) +
  geom_smooth(aes(y = 10^N_pred, x = 10^N), method = "lm", col = "black", lwd = .3) +
  geom_abline(intercept = 0, slope = 1, lwd = .1) +
  facet_grid(direction~interaction, scales = "free_y") +
  labs(x = expression(N[sim]), 
       y = expression(N[GAM]), 
       col = "Scenario") +
  theme_bw() +
  theme(legend.position = "bottom") 
ggsave("figures/figsupp_GAMvSim.png", width = 11.8, height = 9)

# plot standard error from the GAM predictions
ggplot(filter(df_err, Lag == "0")) +
  geom_jitter(aes(y = N_se, 
                  x = interaction, 
                  col = direction), 
              size = 1, alpha = .2,
              position = position_jitterdodge()) +
  geom_boxplot(aes(y = N_se, 
                  x = interaction, 
                 fill = direction), alpha = .6, outlier.shape = NA) +
  facet_wrap(~Process_error, dir = "v", labeller = label_both) +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,3,5)]) +
  scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,3,5)]) +
  labs(x = "", y = "GAM standard error", fill = "Scenario", col = "Scenario") +
  theme(legend.position = "top")

ggline(filter(df_err, Lag == "0"),
       y = "N_se",
       x = "interaction",
       color = "direction",
       facet.by = "Process_error",
       point.size = 2, size = .3,
       add = "mean_sd") +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,3,5)]) +
  facet_wrap(~Process_error, dir = "v", labeller = as_labeller(facet_names)) +
  labs(x = "", y = "GAM standard error", 
       fill = "Scenario", col = "Scenario") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  theme(axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.grid.major.y = element_line(),
        legend.position = "top")

ggsave("figures/figsupp_GAMerror.png", width = 8.5, height = 8.79)

# # the bootstrapped confidence intervals always underestimate the error-propagated intervals
# a <- ggplot(filter(df_err, Lag == "0")) +
#   geom_jitter(aes(y = cilo_se-CI_low, 
#                   x = interaction, 
#                   col = direction), 
#               size = 1, alpha = .2,
#               position = position_jitterdodge()) +
#   geom_boxplot(aes(y = cilo_se-CI_low, 
#                    x = interaction, 
#                    fill = direction), alpha = .6, outlier.shape = NA) +
#   facet_wrap(~Process_error, dir = "v") +
#   labs(title = "Lower confidence interval",
#        x = "", 
#        y = expression(Delta~CI[error]-CI[rlpi]), 
#        fill = "Scenario", col = "Scenario") +
#   theme(legend.position = "bottom")
# 
# b <- ggplot(filter(df_err, Lag == "0")) +
#   geom_jitter(aes(y = cihi_se-CI_high, 
#                   x = interaction, 
#                   col = direction), 
#               size = 1, alpha = .2,
#               position = position_jitterdodge()) +
#   geom_boxplot(aes(y = cihi_se-CI_high, 
#                    x = interaction, 
#                    fill = direction), alpha = .6, outlier.shape = NA) +
#   facet_wrap(~Process_error, dir = "v") +
#   labs(title = "Upper confidence interval",
#        x = "", 
#        y = expression(Delta~CI[error]-CI[rlpi]), 
#        fill = "Scenario", col = "Scenario") +
#   theme(legend.position = "bottom")
# a + b
# ggsave("figures/figsupp_GAMvsBoot_CI.png", width = 17.3, height = 7)