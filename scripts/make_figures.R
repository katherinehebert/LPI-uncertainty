# Script to make figures for the manuscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

## DATA ----

# load results
df0 <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/all_results.RDS")

# set factor levels for plotting
df0$Process_error <- factor(df0$Process_error, levels = c("0", "0.1", "0.2"))
df0$interaction <- gsub("-0.2", "Strong Synchrony", df0$interaction) 
df0$interaction <- gsub("-0.1", "Weak Synchrony", df0$interaction) 
df0$interaction <- gsub("0.1", "Weak Asynchrony", df0$interaction) 
df0$interaction <- gsub("0.2", "Strong Asynchrony", df0$interaction) 
df0$interaction <- gsub("0", "No Synchrony", df0$interaction) 
df0$interaction <- factor(df0$interaction, 
                         levels = c("Strong Asynchrony", "Weak Asynchrony", "No Synchrony", "Weak Synchrony", "Strong Synchrony"))

## PALETTES ----

# colour palette for interactions
colours <- c("Strong Asynchrony" = "#e66101", 
             "Weak Asynchrony" = "#fdb863", 
             "Weak Synchrony"= "#b2abd2", 
             "Strong Synchrony" = "#5e3c99")

# colour palette for process error plots
pal_lags <- rev(c("#22A884FF", "#2A788EFF", "#414487FF"))
pal_error <- c("#08306b", "#2171b5", "#6baed6")
pal_compare <- c("raw" = "blue", "predicted" = "red")

## THEME ----

theme_set(ggpubr::theme_pubr())
format_lpiplots <- list(
  theme(legend.position = "none"),
  labs(x = "", y = "LPI"),
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 2)),
  ylim(c(0.6, 1.6)))

## FIG 1: Overview of the scenario trends

# 1A-C, 2A-C, 3A-C
df <- dplyr::filter(df0, Process_error == "0" & Lag == "0")

# simulated pops
scenarios <- lapply(unique(df$scenario),
                    function(x) readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/simulations/",x,"_l.RDS"))) 
names(scenarios) <- unique(df$scenario)
scenarios <- bind_rows(scenarios, .id = "scenario")
scenarios <- left_join(scenarios, 
                       distinct(subset(df, select = c(scenario, direction, interaction))), by = "scenario")
# simulations
ONE_A <- ggplot(scenarios, 
       aes(x = time, col = direction, group = interaction(popID, scenario))) +
  geom_line(aes(y = N), lwd = .2) +
  labs(x = "", y = "Abundance (N)") + 
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 2)) +
  facet_wrap(~interaction, nrow = 5) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# lpi
ONE_B <- ggplot(df, 
       aes(x = time, col = direction, group = scenario)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, 
                    fill = direction), alpha = .3, lwd = 0) +
    geom_line(aes(y = LPI_final_true), lty = 2, lwd = .2) +
    geom_line(aes(y = LPI_final)) +
    format_lpiplots +
  labs(col = "Direction of change", fill = "Direction of change") +
  facet_wrap(~interaction, nrow = 5) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# put together and save
ONE_A + ONE_B + plot_annotation(tag_levels = 'a')
ggsave("figures/fig2_trendoverview.png", height = 10.7, width = 5.6)


## FIG 2: Accuracy of the trends

# 2a: accuracy
accuracy_plot <- function(df_subset, comparison_variable, colour_variable = "direction") {
  ggplot(df_subset,
         aes(x = get(comparison_variable),
             y = accuracy, 
             group = scenario#,
             #col = get(colour_variable)
             )) +
    geom_hline(yintercept = 0, lwd = 0.3, lty = 2) +
    geom_violin(aes(fill = get(colour_variable),
                    y = accuracy),
                alpha = .5, lwd = .2,
                position = position_dodge(width = .5)) +
    stat_summary(fun.data = mean_se,
                 fun.args = list(mult=1),
                 geom = "pointrange",
                 position = position_dodge(width = .5)) +
    labs(x = "", 
         y = "Bias in the LPI", 
         col = "Direction of change", 
         fill = "Direction of change") +
    theme(legend.position = "top") +
    coord_cartesian(ylim = c(-0.05, 0.05)) 
}
(TWO_A <- accuracy_plot(df, comparison_variable = "interaction"))

# 2b: Percentile of the true LPI within the CI
percentile_plot <- function(df_subset, comparison_variable, colour_variable = "direction") {
  ggplot(data = df_subset, aes(x = get(comparison_variable), 
                               group = scenario
                               )) +
    geom_boxplot(aes(fill = get(colour_variable),
                    y = percentile),
                alpha = .5, lwd = .5, width = .4,
                position = position_dodge(width = .5), outlier.alpha = 0) +
    geom_hline(aes(yintercept = 0.025), lty = 4) +
    geom_hline(aes(yintercept = 0.5), lty = 2) +
    geom_hline(aes(yintercept = 0.975), lty = 4) +
    labs(y = expression(mu~Percentile), 
         x = "", 
         col = "Direction of change", 
         fill = "Direction of change") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(-0.1,1.1)) +
    scale_y_continuous(breaks = c(0.025, 0.5, 0.975))
}
(TWO_B <- percentile_plot(df, comparison_variable = "interaction"))

# put together and save
TWO_A / TWO_B + plot_annotation(tag_levels = 'a')
ggsave("figures/fig3_accuracy.png", width = 10.7, height = 5.6)



## FIG 3: Precision of the trends

#3a: Does the residual error match the gam error??
THREE_A <- ggplot(df, aes(x = interaction, y = residual_error_sd, fill = direction, group = scenario)) +
  geom_violin(alpha = .5, lwd = .2, width = .5,
              position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0.05, lty = 2) +
  labs(x = "", 
       y = "Standard deviation of the\nresidual error distribution",
       col = "Direction of change", 
       fill = "Direction of change") +
  stat_summary(aes(y = residual_error_sd),
               fun.data = mean_se,
               fun.args = list(mult=1),
               geom = "pointrange",
               position = position_dodge(width = .5)) +
  coord_cartesian(ylim = c(0, 0.06))

#3b: Overlap between raw and smoothed(?) distributions

overlap_plot <- function(df_subset, comparison_variable, colour_variable = "direction") {
  ggplot(data = df_subset, aes(x = get(comparison_variable), 
                               group = scenario#,
                               #col = get(colour_variable)
                               )) +
    geom_violin(aes(fill = get(colour_variable),
                    y = 100*overlap),
                alpha = .5, lwd = .2, width = .65,
                position = position_dodge(width = .5)) +
    stat_summary(aes(y = 100*overlap),
                 fun.data = mean_se,
                 fun.args = list(mult=1),
                 geom = "pointrange",
                 position = position_dodge(width = .5)) +
    labs(y = "Overlap (%)", 
         x = "",  
         col = "Direction of change", 
         fill = "Direction of change"
         ) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0,50))
}
THREE_B <- overlap_plot(df, comparison_variable = "interaction")

THREE_A / THREE_B + plot_annotation(tag_levels = "a")
ggsave("figures/fig4_precision.png", width = 10.7, height = 5.6)


## SUPPLEMENTARY

# reasonable precision plots
# process error
# - lag should go in supplementary

# LAG #######

# load results
df <- dplyr::filter(df0, Process_error == "0")

# 4a: accuracy
(FOUR_A <- accuracy_plot(filter(df, interaction != "No Synchrony"), 
                         colour_variable = "Lag", 
                         comparison_variable = "direction") + 
    facet_wrap(~interaction, ncol = 4) + 
    scale_fill_manual(values = pal_lags) +
  labs(fill = "Lag", color = "Lag"))
# 4b: percentile
(FOUR_B <- percentile_plot(filter(df, interaction != "No Synchrony"), 
                           colour_variable = "Lag", 
                           comparison_variable = "direction") + 
    facet_wrap(~interaction, ncol = 4) + 
    scale_fill_manual(values = pal_lags) +
  labs(fill = "Lag"))

#5a: Does the residual error match the gam error??
(FOUR_C <- ggplot(filter(df, interaction != "No Synchrony"), 
                  aes(x = direction, 
                      y = residual_error_sd, 
                      fill = Lag, 
                      group = scenario)) +
  geom_violin(alpha = .5, lwd = .2, width = .5,
              position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0.05, lty = 2) + 
  facet_wrap(~interaction, ncol = 4) + 
  scale_fill_manual(values = pal_lags) +
  labs(x = "", 
       y = "Standard deviation of the\nresidual error distribution",
       col = "Lag", 
       fill = "Lag") +
  stat_summary(aes(y = residual_error_sd),
               fun.data = mean_se,
               fun.args = list(mult=1),
               geom = "pointrange",
               position = position_dodge(width = .5)) +
  coord_cartesian(ylim = c(0, 0.06)) +
    theme(legend.position = "none"))
#5b: Overlap between raw and smoothed(?) distributions
(FOUR_D <- overlap_plot(filter(df, interaction != "No Synchrony"), 
                       comparison_variable = "direction",
                       colour_variable = "Lag") +
  facet_wrap(~interaction, ncol = 4) + 
  scale_fill_manual(values = pal_lags))


FOUR_A / FOUR_B / FOUR_C / FOUR_D + plot_annotation(tag_levels = "a")
ggsave("figures/fig5_lag.png", width = 10.7, height = 12)


### Supplementary figures ----

#### Process error 

df <- dplyr::filter(df0, Lag == "0")

# S2a: accuracy
(S2_A <- accuracy_plot(df, 
                       colour_variable = "Process_error", 
                         comparison_variable = "direction") + 
    facet_wrap(~interaction, ncol = 5) + 
    scale_fill_brewer() +
    labs(fill = "Process error", color = "Process error"))
# S2b: percentile
(S2_B <- percentile_plot(df, 
                         colour_variable = "Process_error", 
                           comparison_variable = "direction") + 
    facet_wrap(~interaction, ncol = 5) + 
    scale_fill_brewer() +
    labs(fill = "Process error"))

# S2c: Does the residual error match the gam error??
(S2_C <- ggplot(df,
                  aes(x = direction, 
                      y = residual_error_sd, 
                      fill = Process_error, 
                      group = scenario)) +
    geom_violin(alpha = .5, lwd = .2, width = .5,
                position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0.05, lty = 2) + 
    facet_wrap(~interaction, ncol = 5) + 
    scale_fill_brewer() +
    labs(x = "", 
         y = "Standard deviation of the\nresidual error distribution",
         col = "Process error", 
         fill = "Process error") +
    stat_summary(aes(y = residual_error_sd),
                 fun.data = mean_se,
                 fun.args = list(mult=1),
                 geom = "pointrange",
                 position = position_dodge(width = .5)) +
    coord_cartesian(ylim = c(0, 0.06)) +
    theme(legend.position = "none"))
#S2d: Overlap between raw and smoothed(?) distributions
(S2_D <- overlap_plot(df, 
                        comparison_variable = "direction",
                        colour_variable = "Process_error") +
    facet_wrap(~interaction, ncol = 5) + 
    scale_fill_brewer()) 

S2_A / S2_B / S2_C / S2_D + plot_annotation(tag_levels = "a")
ggsave("figures/figS2_processerror.png", width = 10.7, height = 12)


#### Process error 

df <- df0
df$Lag <- paste0("Lag-", df$Lag)

# S3: accuracy
(S3 <- accuracy_plot(filter(df, interaction != "No Synchrony"), 
                     colour_variable = "Process_error", 
                     comparison_variable = "Lag") + 
    facet_grid(interaction~direction) + 
    scale_fill_brewer() +
    labs(fill = "Process error", color = "Process error") +
    coord_cartesian(ylim = c(-0.04, 0.04))) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line())
ggsave("figures/figS3_processerror_lag.png", width = 12.3, height = 10)


# S4: percentile
(S4 <- percentile_plot(filter(df, interaction != "No Synchrony"), 
                       colour_variable = "Process_error", 
                       comparison_variable = "Lag") + 
    facet_grid(interaction~direction) + 
    scale_fill_brewer() +
    labs(fill = "Process error") +
    theme(legend.position = "top"))
ggsave("figures/figS4_processerror_lag.png", width = 12.3, height = 10)


# S5: Does the residual error match the gam error??
(S5 <- ggplot(filter(df, interaction != "No Synchrony"),
              aes(x = Lag, 
                  y = residual_error_sd, 
                  fill = Process_error, 
                  group = scenario)) +
    geom_violin(alpha = .5, lwd = .2, width = .5,
                position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0.05, lty = 2) + 
    facet_grid(interaction~direction) + 
    scale_fill_brewer() +
    labs(x = "", 
         y = "Standard deviation of the\nresidual error distribution",
         col = "Process error", 
         fill = "Process error") +
    stat_summary(aes(y = residual_error_sd),
                 fun.data = mean_se,
                 fun.args = list(mult=1),
                 geom = "pointrange",
                 position = position_dodge(width = .5)) +
    coord_cartesian(ylim = c(0, 0.06)) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line()))
ggsave("figures/figS5_processerror_lag.png", width = 12.3, height = 10)

# S6: Overlap between raw and smoothed(?) distributions
(S2_D <- overlap_plot(filter(df, interaction != "No Synchrony"),
                      comparison_variable = "Lag",
                      colour_variable = "Process_error") +
    facet_grid(interaction~direction) + 
    scale_fill_brewer() +
    theme(panel.grid.major.y = element_line())) 
ggsave("figures/figS6_processerror_lag.png", width = 12.3, height = 10)
