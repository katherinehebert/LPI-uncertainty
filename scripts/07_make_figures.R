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
df0 <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/all_results.RDS")
# df2 <- readRDS("outputs/all_uncertaintypropagation.RDS")
# df0 <- left_join(df1, df2, by = c("scenario", "time"))

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

facet_names <- c(
  `0` = "Process ε = 0",
  `0.1` = "Process ε = 0.1",
  `0.2` = "Process ε = 0.2",
  `Strong Asynchrony` = "Strong Asynchrony", 
  `Weak Asynchrony` = "Weak Asynchrony", 
  `No Synchrony` = "No Synchrony",
  `Weak Synchrony`= "Weak Synchrony", 
  `Strong Synchrony` = "Strong Synchrony"
)

# colour palette for process error plots
pal_lags <- rev(c("#22A884FF", "#2A788EFF", "#414487FF"))
pal_error <- c("#08306b", "#2171b5", "#6baed6")
pal_compare <- c("raw" = "blue", "predicted" = "red")

## THEME ----

theme_set(ggpubr::theme_pubr())
format_lpiplots <- list(
  theme(legend.position = "none"),
  labs(x = "", y = "LPI"),
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)),
  ylim(c(0.3, 2.3)))

## FIG 1: Overview of the scenario trends ----

# 1A-C, 2A-C, 3A-C
df <- dplyr::filter(df0, Process_error == "0.1" & Lag == "0")
df$lpi_variance[which(df$time == 1)] = 0

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
                aes(x = time-1, col = direction, group = interaction(popID, scenario))) +
  geom_line(aes(y = N), lwd = .3) +
  labs(x = "Time", y = "Abundance (N)", col = "Trend") + 
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
  facet_wrap(~interaction, ncol = 5) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.grid.major.y = element_line(),
        panel.spacing.x = unit(4, "mm")) +
  coord_cartesian(ylim = c(0, 2750))

(ONE_B <- ggplot(df, 
                 aes(x = time-1, col = direction, group = scenario)) +
    geom_ribbon(aes(ymin = CI_low, 
                    ymax = CI_high, 
                    fill = direction), 
                alpha = .9, lwd = 0) +
    geom_ribbon(aes(ymin = lpi_correction - 1.96*sqrt(lpi_variance), 
                    ymax = lpi_correction + 1.96*sqrt(lpi_variance), 
                    col = direction, 
                    fill = direction), 
                alpha = .3, lwd = .1) +
    geom_line(aes(y = LPI_final, lty = "Original"), col = "black", lwd = .4) +
    geom_line(aes(y = lpi_correction, lty = "Corrected"), col = "black", lwd = .4) +
    format_lpiplots +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    scale_fill_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    scale_linetype_manual(values = c("Original" = 2, "Corrected" = 1)) +
    labs(x = "Time", col = "Trend", fill = "Trend", linetype = "LPI") +
    facet_wrap(~interaction, ncol = 5) +
    theme(legend.position = "top",
          axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12),
          panel.grid.major.y = element_line(),
          panel.spacing.x = unit(4, "mm"))
)

# put together and save
(ONE_A / ONE_B + plot_annotation(tag_levels = 'a'))
ggsave("figures/fig1_trendoverview.png", width = 9, height = 8)

#### remove time = 1 ----
# bc it is just the baseline

df0 <- dplyr::filter(df0, time != 1)

## FIG 2 & 3 : Accuracy of the trend and the uncertainty interval ----

df <- dplyr::filter(df0, Lag == 0)

(FIG2_A <- ggline(df,
                  "Process_error",
                  "percentile_LPIcorrected",
                  facet.by = "interaction",
                  color = "direction",
                  point.size = 2, size = .5, 
                  add = c("mean_sd")) +
    facet_wrap(~interaction, ncol = 5, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    geom_hline(aes(yintercept = 0.025), lty = 4, alpha = .4) +
    geom_hline(aes(yintercept = 0.5), lty = 2, alpha = .4) +
    geom_hline(aes(yintercept = 0.975), lty = 4, alpha = .4) +
    labs(y = "Percentile of the corrected LPI \nin the original CI", #expression(mu~Percentile),
         x = "",
         col = "Trend",
         fill = "Trend") +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12),
          legend.position = "none") +
    coord_cartesian(ylim = c(-0.1,1.1)) +
    scale_y_continuous(breaks = c(0.025, 0.5, 0.975)))

(FIG2_B <- ggline(df,
                  "Process_error",
                  "CI_width_reldiff",
                  facet.by = "interaction",
                  color = "direction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, ncol = 5, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = " ",
         y = "Unrepresented uncertainty") +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12),
          panel.grid.major.y = element_line(),
          legend.position = "none"))

(FIG3_A <- ggline(df, 
                  "Process_error", 
                  "lpi_bias",
                  color = "direction", 
                  point.size = 2, size = .5,
                  facet.by = "interaction",
                  add = c("mean_sd")) +
    facet_wrap(~interaction, dir = "h", ncol = 5) + 
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = " ", 
         y = "Uncertainty bias of the LPI") +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12),
          panel.grid.minor.y = element_line(linewidth = .2),
          panel.grid.major.y = element_line(),
          legend.position = "top",
          panel.spacing.x = unit(4, "mm")) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG3_B <- ggline(df, 
                  "Process_error", 
                  "lpi_variance",
                  color = "direction", 
                  point.size = 2, size = .5,
                  facet.by = "interaction",
                  add = c("mean_sd")) +
    facet_wrap(~interaction, dir = "h", nrow = 1) + #+
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = " ", 
         y = "Variance of the LPI") +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 12),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(linewidth = .2),
          legend.position = "top",
          panel.spacing.x = unit(4, "mm")))

(FIG3_A / FIG2_A + plot_annotation(tag_levels = "a"))
ggsave("figures/fig2_accuracy.png", width = 9, height = 8)

(FIG3_B / FIG2_B + plot_annotation(tag_levels = "a")) 
ggsave("figures/fig3_uncertainty.png", width = 9, height = 8)

## FIG 4: Demonstration of uncertainty of population trends vs. pop size ----

# equation 2 ----
# expectation of the growth rate trend
eq2 <- function(N, sigma_measure = 0){
  d <- c()
  for(t in 2:length(N)){
    d[t] <- log10(N[t]/N[t-1]) + (sigma_measure^2)/(2*(N[t-1]^2 - N[t]^2))
  }
  return(d)
}

# equation 2
# just the growth rate part
eq2_donly <- function(N){
  eq2dt <- c()
  for(t in 2:length(N)){
    eq2dt[t] <- log10(N[t]/N[t-1])
  }
  return(eq2dt)
}

eq2_varonly <- function(N, sigma_measure = 0){
  eq2var <- c()
  for(t in 2:length(N)){
    eq2var[t] <- (sigma_measure^2)/(2*((N[t-1])^2 - (N[t])^2))
  }
  return(eq2var)
}

N = 2:200
N2 = N-1
NN = cbind(N, N2)
correction <- apply(NN, 1, eq2_varonly, sigma_measure = 5)[2,]
NNC = data.frame(
  "N0" = N,
  "N1" = N2,
  "correction" = correction,
  "Trend" = "Decline ( N - 1 )"
)

Nd2 = N+1
NNd = cbind(N, Nd2)
correction <- apply(NNd, 1, eq2_varonly, sigma_measure = 5)[2,]
NNCd = data.frame(
  "N0" = N,
  "N1" = Nd2,
  "correction" = correction,
  "Trend" = "Growth ( N + 1 )"
)
NNC = rbind(NNC, NNCd) %>% as.data.frame()


ggplot(data = NNC) +
  geom_line(aes(x = N0, y = correction, col = Trend), lwd = .8) +
  scale_color_manual(values = pal_locuszoom("default")(6)[c(1,3)]) +
  geom_hline(yintercept = 0, lwd = .2, lty = 2) +
  labs(y = expression(Measurement~uncertainty~correction~to~"E["~d[t]~"]"), 
       x = "Population abundance (N)") +
  theme_pubr()
ggsave("figures/fig4_demo_eq2.png", width = 5.5, height = 4.5)


## SUPPLEMENTARY FIGURES ####

# FIG S: LAG ----

# load results
df <- dplyr::filter(df0, Process_error == "0" & interaction != "No Synchrony")

(FIG4_A <- ggline(df, 
                  "Lag", 
                  "CI_width_reldiff",
                  color = "direction",
                  facet.by = "interaction",
                  point.size = 2, size = .5,
                  add = c("mean_sd")) +
    facet_wrap(~interaction, nrow =  1) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Unrepresented uncertainty") +
    scale_x_discrete(labels = function(x) paste0("Lag-", x)) +
    theme(axis.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          panel.grid.major.y = element_line())+
    #coord_cartesian(ylim = c(-0.05, 0.05)) +
    geom_hline(yintercept = 0, lwd = .2, lty = 2))
(FIG4_B <- ggline(df, 
                  "Lag", 
                  "percentile_LPIcorrected",
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
    labs(y = "Percentile of the\ncorrected LPI", 
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
(FIG5_C / FIG4_B + plot_annotation(tag_levels = "a")) 
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
    coord_cartesian(ylim = c(-0.07, 0.01)) +
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
(FIG5_D / FIG4_A + plot_annotation(tag_levels = "a")) 
ggsave("figures/fig5_lag_uncertainty.png", width = 8.56, height = 6.77)

# FIG S: LAG + ERROR ----

df <- dplyr::filter(df0, interaction != "No Synchrony")

(FIGSX_A <- ggline(df, 
                   "Lag", 
                   "CI_width_diff",
                   color = "direction",
                   facet.by = c("interaction", "Process_error"),
                   point.size = 2, size = .3,
                   add = c("mean_sd")) +
    facet_wrap(~Process_error+interaction, nrow = 3, ncol = 4, labeller = as_labeller(facet_names)) +
    scale_color_manual(values = pal_locuszoom("default")(6)[c(1,5,3)]) +
    labs(color = "Trend",
         x = "Covariance lag", 
         y = "Unrepresented uncertainty") +
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
                   "percentile_LPIcorrected",
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
    labs(y = "Percentile of the\ncorrected LPI", #expression(mu~Percentile), 
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
    coord_cartesian(ylim = c(-0.5, 0.1)) +
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