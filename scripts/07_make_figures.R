# Script to make figures for the manuscript

# LOAD LIBRARIES ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(ggsci)

## PREPARE DATA ----

# load results
df1 <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/all_results.RDS")
df2 <- readRDS("outputs/all_uncertaintypropagation.RDS")
df0 <- left_join(df1, df2, by = c("scenario", "time"))

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
#   scale_color_viridis_d(begin = .3, end = .7, option = "A") 

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
  ylim(c(0.5, 2)))

## FIG 1: Overview of the scenario trends ----

# 1A-C, 2A-C, 3A-C
df <- dplyr::filter(df0, Process_error == "0" & Lag == "0")

# load simulated populations
scenarios <- lapply(unique(df$scenario),
                    function(x) readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/",x,"_l_dtchain.RDS"))) 
names(scenarios) <- unique(df$scenario)
scenarios <- bind_rows(scenarios, .id = "scenario")
scenarios <- left_join(scenarios, 
                       distinct(subset(df, select = c(scenario, direction, interaction))), by = "scenario")
scenarios$dt_chain[which(scenarios$time == 1)] <- 0

# plot the simulated population trends
ONE_A <- ggplot(scenarios, 
       aes(x = time, col = direction, group = interaction(popID, scenario))) +
  geom_line(aes(y = N), lwd = .2) +
  labs(x = "Time", y = "Abundance (N)", col = "Trend") + 
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 2)) +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  facet_wrap(~interaction, nrow = 5) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
# plot the corresponding lpi trends
ONE_B <- ggplot(df, 
       aes(x = time, col = direction, group = scenario)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high, 
                    fill = direction), alpha = .3, lwd = 0) +
    geom_line(aes(y = LPI_final_true), lty = 2, lwd = .2) +
    geom_line(aes(y = LPI_final)) +
    format_lpiplots +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  labs(x = "Time", col = "Trend", fill = "Trend") +
  facet_wrap(~interaction, nrow = 5) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

ONE_C <- ggplot(df, 
                aes(x = time, col = direction, group = scenario)) +
  geom_line(aes(y = lpi_correction), lty = 2, lwd = .2) +
  geom_line(aes(y = lpi_nocorrection)) +
  format_lpiplots +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  labs(col = "Trend", fill = "Trend") +
  facet_wrap(~interaction, ncol = 5) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# put together and save
(ONE_A + ONE_B + plot_annotation(tag_levels = 'a') )
ggsave("figures/fig1_trendoverview.png", height = 7, width = 6)

# plot the simulated population trends
ggplot(scenarios, 
                aes(x = time, col = direction, group = interaction(popID, scenario))) +
  geom_line(aes(y = N), lwd = .2) +
  labs(x = "", y = "Abundance (N)") + 
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 2)) +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  facet_wrap(~interaction, ncol = 5) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
ggsave("figures/presentation_trendoverview.png", width = 8.89, height = 2.2)

# ggplot(dplyr::filter(scenarios, scenario == "scenario1C"),
#        aes(x = time, col = popID)) +
#   geom_line(aes(y = N), lwd = .2) +
#   ggpubr::theme_transparent() +
#   theme(legend.position = "none")
# ggsave("figures/presentation_titleslide.png", width = 13.4, height = 5.45)


#### remove time = 1 ----
# bc it is just the baseline

df0 <- dplyr::filter(df0, time != 1)


## FIG 2: Accuracy of the trends ----

df <- dplyr::filter(df0, Lag == 0, Process_error == 0)

(FIG2_A <- ggline(df, 
                  "interaction", 
                  "accuracy",
                  color = "direction", 
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = " ", 
         y = "Bias from the expected LPI") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line(),
          legend.position = "none") +
    coord_cartesian(ylim = c(-0.05, 0.05)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG2_B <- ggline(df, 
                  "interaction", 
                  "percentile",
                  color = "direction", 
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    geom_hline(aes(yintercept = 0.025), lty = 4, alpha = .4) +
    geom_hline(aes(yintercept = 0.5), lty = 2, alpha = .4) +
    geom_hline(aes(yintercept = 0.975), lty = 4, alpha = .4) +
    labs(y = "Percentile of the expected LPI", #expression(mu~Percentile), 
         x = "", 
         col = "Trend", 
         fill = "Trend") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = "bottom") +
    coord_cartesian(ylim = c(-0.1,1.1)) +
    scale_y_continuous(breaks = c(0.025, 0.5, 0.975)))
FIG2_leg <- ggpubr::get_legend(FIG2_B) %>% as_ggplot()
((FIG2_leg / (FIG2_A + (FIG2_B + theme(legend.position = "none")))) + 
    plot_annotation(tag_levels = "a")) +
  plot_layout(height = c(1,5))
ggsave("figures/fig2_accuracy.png", width = 12.4, height = 5.14)


## FIG 3: Uncertainty (propagated) ----

df <- dplyr::filter(df0, Lag == 0)

facet_names <- c(
  `0` = "Process ε = 0",
  `0.1` = "Process ε = 0.1",
  `0.2` = "Process ε = 0.2",
  `Strong Asynchrony` = "Strong Asynchrony", 
  `Weak Asynchrony` = "Weak Asynchrony", 
  `Weak Synchrony`= "Weak Synchrony", 
  `Strong Synchrony` = "Strong Synchrony"
)

(FIG3_A <- ggline(df, 
               "interaction", 
               "lpi_bias",
               color = "direction", 
               point.size = 2, size = .5,
               facet.by = "Process_error",
               add = c("mean_sd")) +
  facet_wrap(~Process_error, dir = "h", labeller = as_labeller(facet_names)) + 
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  labs(color = "Trend",
       x = " ", 
       y = "Uncertainty bias of the LPI") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10),
                   breaks = function(x){x[c(TRUE, FALSE)]}) +
  theme(axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.grid.major.y = element_line(),
        legend.position = "top",
        panel.spacing.x = unit(4, "mm")) +
  #coord_cartesian(ylim = c(-0.3, 0.1)) +
  geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG3_B <- ggline(df, 
                "interaction", 
                "lpi_variance",
                color = "direction", 
                point.size = 2, size = .5,
                facet.by = "Process_error",
                add = c("mean_sd")) +
    facet_wrap(~Process_error, dir = "h", labeller = as_labeller(facet_names)) + #+
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = " ", 
         y = "Variance of the LPI") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10),
                     breaks = function(x){x[c(TRUE, FALSE)]}) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line(),
          legend.position = "none",
          panel.spacing.x = unit(4, "mm")) +
    coord_cartesian(ylim = c(0, 0.5)))
(FIG3_A / FIG3_B + plot_annotation(tag_levels = "a")) 
ggsave("figures/fig3_uncertainty.png", width = 10.9, height = 7.7)


## SUPPLEMENTARY FIGURES ####

# FIG S: LAG ----

# load results
df <- dplyr::filter(df0, Process_error == "0" & interaction != "No Synchrony")

(FIG4_A <- ggline(df, 
                  "Lag", 
                  "accuracy",
                  color = "direction",
                  facet.by = "interaction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, nrow =  1) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Bias from the\nexpected LPI") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line())+
    #coord_cartesian(ylim = c(-0.05, 0.05)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG4_B <- ggline(df, 
                  "Lag", 
                  "percentile",
                  color = "direction", 
                  facet.by = "interaction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, nrow =  1) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    geom_hline(aes(yintercept = 0.025), lty = 4, alpha = .4) +
    geom_hline(aes(yintercept = 0.5), lty = 2, alpha = .4) +
    geom_hline(aes(yintercept = 0.975), lty = 4, alpha = .4) +
    labs(y = "Percentile of the\nexpected LPI", 
         x = "Covariance lag", 
         col = "Trend", 
         fill = "Trend") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = "none") +
    #coord_cartesian(ylim = c(-0.1,1.1)) +
    scale_y_continuous(breaks = c(0.025, 0.5, 0.975)))
(FIG4_A / FIG4_B + plot_annotation(tag_levels = "a")) 
ggsave("figures/fig4_lag_accuracy.png", width = 8.56, height = 6.77)

(FIG5_C <- ggline(df, 
                  "Lag", 
                  "lpi_bias",
                  color = "direction",
                  facet.by = "interaction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, nrow =  1) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Uncertainty bias of the LPI") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line())+#,
    #legend.position = "none") +
    coord_cartesian(ylim = c(-0.6, 0.1)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG5_D <- ggline(df, 
                  "Lag", 
                  "lpi_variance",
                  color = "direction", 
                  facet.by = "interaction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, nrow =  1) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(y = "Variance of the LPI", #expression(mu~Percentile), 
         x = "Covariance lag", 
         col = "Direction\n of change", 
         fill = "Trend") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = "none") #+
    #coord_cartesian(ylim = c(0,0.1))
  ) 
(FIG5_C / FIG5_D + plot_annotation(tag_levels = "a")) 
ggsave("figures/fig5_lag_uncertainty.png", width = 8.56, height = 6.77)

# FIG S: LAG + ERROR ----

df <- dplyr::filter(df0, interaction != "No Synchrony")

(FIGSX_A <- ggline(df, 
                  "Lag", 
                  "accuracy",
                  color = "direction",
                  facet.by = c("interaction", "Process_error"),
                  point.size = 2, size = .3,
                  add = c("mean_sd")) +
    facet_wrap(~Process_error+interaction, nrow = 3, ncol = 4, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Bias from the\nexpected LPI") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line())+
    #coord_cartesian(ylim = c(-0.05, 0.05)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
ggsave("figures/figSX_lag_accuracy.png", width = 8.56, height = 8)

(FIGSX_B <- ggline(df, 
                  "Lag", 
                  "percentile",
                  color = "direction", 
                  facet.by = c("interaction", "Process_error"),
                  point.size = 2, size = .3,
                  add = c("mean_sd")) +
    facet_wrap(~Process_error+interaction, nrow = 3, ncol = 4, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    geom_hline(aes(yintercept = 0.025), lty = 4, alpha = .4) +
    geom_hline(aes(yintercept = 0.5), lty = 2, alpha = .4) +
    geom_hline(aes(yintercept = 0.975), lty = 4, alpha = .4) +
    labs(y = "Percentile of the\nexpected LPI", #expression(mu~Percentile), 
         x = "Covariance lag", 
         col = "Trend", 
         fill = "Trend") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = "top") +
    coord_cartesian(ylim = c(-0.1,1.1)) +
    scale_y_continuous(breaks = c(0.025, 0.5, 0.975)))
ggsave("figures/figSX_lag_percentile.png", width = 8.56, height = 8)

(FIGSX_C <- ggline(df, 
                  "Lag", 
                  "lpi_bias",
                  color = "direction",
                  facet.by = c("interaction", "Process_error"),
                  point.size = 2, size = .3,
                  add = c("mean_sd")) +
    facet_wrap(~Process_error+interaction, nrow = 3, ncol = 4, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Uncertainty bias of the LPI") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line())+
    coord_cartesian(ylim = c(-0.6, 0.1)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
ggsave("figures/figSX_lag_uncertaintybias.png", width = 8.56, height = 8)

(FIGSX_D <- ggline(df, 
                  "Lag", 
                  "lpi_variance",
                  color = "direction", 
                  facet.by = c("interaction", "Process_error"),
                  point.size = 2, size = .3,
                  add = c("mean_sd")) +
    facet_wrap(~Process_error+interaction, nrow = 3, ncol = 4, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(y = "Variance of the LPI", #expression(mu~Percentile), 
         x = "Covariance lag", 
         col = "Trend", 
         fill = "Trend") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line(),
          legend.position = "top") #+
    #coord_cartesian(ylim = c(0,0.12))
  ) 
ggsave("figures/figSX_lag_variance.png", width = 8.56, height = 8)

# FIG S: GAM error plots ----

(FIGSX_E <- ggline(filter(df0, interaction != "No Synchrony" & Lag == "0"), 
          "interaction", 
          "residual_error_sd",
          color = "direction", 
          facet.by = c("Process_error"),
          point.size = 2, size = .3,
          add = c("mean_sd")) +
    facet_wrap(~Process_error, dir = "v", labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(y = "Standard deviation of\nthe GAM residual error", #expression(mu~Percentile), 
         x = "", 
         col = "Trend") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line(),
          legend.position = "none") +
    coord_cartesian(ylim = c(0,0.2)) +
  geom_hline(yintercept = 0.05, lty = 2))
#, width = 8.5, height = 8.79

df0$overlap100 <- df0$overlap*100
(FIGSX_F <- ggline(filter(df0, interaction != "No Synchrony" & Lag == "0"), 
                   "interaction", 
                   "overlap100",
                   color = "direction", 
                   facet.by = c("Process_error"),
                   point.size = 2, size = .3,
                   add = c("mean_sd")) +
    facet_wrap(~Process_error, dir = "v", labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    labs(y = "Overlap between raw and smoothed\ngrowth rate distributions (%)", #expression(mu~Percentile), 
         x = "", 
         col = "Trend") +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line(),
          legend.position = "right")) +
  coord_cartesian(ylim = c(0, 20))

FIGSX_E + FIGSX_F + plot_annotation(tag_levels = "a")
ggsave("figures/figsupp_GAMerror.png", width = 10.2, height = 6.98)


# FIG S: GAM predictions vs. simulation plots ----

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
#results$scenario <- uncertainty$scenario
df_err <- inner_join(results, df)

## comparing GAM predictions to simulated population sizes

# plot difference between prediction and simulation
ggplot(filter(df_err, Lag == "0" & time != 1)) +
  geom_jitter(aes(y = (10^N_pred)-(10^N), 
                  x = interaction, 
                  col = direction), size = 1, alpha = .2, position = position_jitterdodge()) +
  facet_wrap(~ Process_error, dir = "v", labeller = as_labeller(facet_names)) +
  labs(x = "", y = expression(N[sim]~Delta~N[GAM]), col = "Trend") +
  theme(legend.position = "top") +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  ggpubr::theme_pubr() +
  geom_hline(yintercept = 0, lwd = .3, lty = 2)
ggsave("figures/figsupp_GAMpredictions.png", width = 6, height = 6)

# plot standard error from the GAM predictions
ggline(filter(df_err, Lag == "0"),
       y = "N_se",
       x = "interaction",
       color = "direction",
       facet.by = "Process_error",
       point.size = 2, size = .3,
       add = "mean_sd") +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  facet_wrap(~Process_error, dir = "v", labeller = as_labeller(facet_names)) +
  labs(x = "", y = "GAM standard error", 
       fill = "Scenario", col = "Scenario") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  theme(axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.grid.major.y = element_line(),
        legend.position = "top")
#ggsave("figures/figsupp_.png", width = 8.5, height = 8.79)