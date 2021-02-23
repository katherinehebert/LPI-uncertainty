# Script to compare results of the simulated population dynamics scenarios

# load packages and prep environment
library(ggplot2)
library(dplyr)
theme_set(ggpubr::theme_pubr())

# data prep for plotting ----

# import results of each scenario's LPI
lpi <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_lpi.RDS")), readRDS)
names(lpi) <- gsub("_lpi.RDS", "", list.files(path = "outputs/", pattern = "_lpi.RDS"))
# and bind into one data frame
lpi <- dplyr::bind_rows(lpi, .id = "scenario")

# import parameter tables for each scenario
params <- lapply(paste0("simulations/", list.files(path = "simulations/", pattern = "_params.RDS")), readRDS)
params <- lapply(params, data.table::transpose) %>%
  lapply(function(x){
    colnames(x) <- gsub(" ", "_", x[1,])
    x <- x[-c(1,3),]
    })
names(params) <- gsub("_params.RDS", "", list.files(path = "simulations/", pattern = "_params.RDS"))
# and bind into one data frame
params <- dplyr::bind_rows(params, .id = "scenario")

# make table for categories from carrying capacity scenarios
K_scenarios <- data.frame("scenario" = params$scenario, "direction" = NA)
K_scenarios$direction[which(K_scenarios$scenario %in% c("scenario01", "scenario04",
                                                        "scenario07", "scenario10",
                                                        "scenario13", "scenario16",
                                                        "scenario19"))]   <- "decline"           
K_scenarios$direction[which(K_scenarios$scenario %in% c("scenario02", "scenario05",
                                                        "scenario08", "scenario11",
                                                        "scenario14", "scenario17",
                                                        "scenario20"))]   <- "stable"   
K_scenarios$direction[which(K_scenarios$scenario %in% c("scenario03", "scenario06",
                                                        "scenario09", "scenario12",
                                                        "scenario15", "scenario18",
                                                        "scenario21"))]   <- "growth" 

# join all tables together
df <- dplyr::left_join(lpi, K_scenarios) %>% dplyr::left_join(params)

# calculate LPI accuracy as % difference [(estimated - true)/true * 100]
df$accuracy_boot <- ((df$LPI_boot - df$LPI_true)/df$LPI_true)*100
# does the true LPI fall within the 95% confidence interval?
if(df$LPI_true < df$cihi_boot && df$LPI_true > df$cilo_boot){
  df$precision_boot <- "yes"} else {
    df$precision_boot <- "no"
  }
# interval width divided by the lpi value
df$interval_width <- (df$cihi_boot - df$cilo_boot)/df$LPI_boot

# format columns for plotting
df$Lag <- factor(df$Lag, levels = c("0", "1", "2"))
df$direction <- factor(df$direction, levels = c("decline", "stable", "growth"))
colnames(df)[11] <- "N0"
colnames(df)[12] <- "lambda"
colnames(df)[13] <- "interaction"


# plotting exploration ----

ggplot(df, aes(x = direction, y = accuracy_boot, col = interaction)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  xlab("") +
  facet_wrap(~ Lag)

# effect of direction only
ggplot(filter(df, Lag == "0" & interaction == 0), 
       aes(y = direction, x = accuracy_boot, col = direction)) +
  geom_jitter(width = .1, alpha = .3) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  labs(y = "", x = "Accuracy (% difference in LPI)") +
  geom_vline(xintercept = 0, lty = 2, lwd = .2) +
  coord_cartesian(xlim = c(-20,20)) +
  theme(legend.position = "none")

# effect of interaction only
ggplot(filter(df, Lag == "0" & direction == "stable"), 
       aes(x = interaction, y = accuracy_boot, col = interaction)) +
  geom_jitter(width = .1, alpha = .3) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  labs(x = "Interaction level", y = "Accuracy (% difference in LPI)") +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2")

# effect of lag and direction only
ggplot(filter(df), 
       aes(y = Lag, x = accuracy_boot, col = Lag)) +
  geom_jitter(width = .1, alpha = .1) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) +
  geom_vline(xintercept = 0, lty = 2, lwd = .2) +
  labs(y = "Lag length (time steps)", x = "Accuracy (% difference in LPI)") +
  coord_cartesian(xlim = c(-30,30)) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~direction, dir = "v")


ggplot(filter(df, Lag == "0" & interaction == 0), 
       aes(group = direction, x = accuracy_boot, fill = direction)) +
  geom_density(stat = "density", alpha = .5, lwd = .1,outline.type = "upper") +
  labs(x = "Density", x = "Accuracy (% difference in LPI)")


# no lag. effect of covariance and direction
ggplot(filter(df, Lag == "0"), aes(y = direction, x = accuracy_boot, col = interaction)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  geom_vline(xintercept = 0, lty = 2, lwd = .2) +
  labs(y = "", x = "Accuracy (% difference in LPI)") +
  coord_cartesian(xlim = c(-12,12)) +
  scale_color_brewer(palette = "Dark2")


ggplot(df, aes(x = direction, y = interval_width, col = interaction)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  xlab("") +
  facet_wrap(~ Lag)

ggplot(filter(df, Lag == "0"), aes(x = direction, y = interval_width, col = interaction)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", position = position_dodge(width = .5)) + 
  xlab("")

ggplot(df, aes(x = Lag, y = interval_width)) +
  geom_point(col = "grey") +
  stat_summary(fun.data = mean_se, fun.args = list(mult=1),
               geom = "pointrange", , position = position_dodge(width = .5)) + 
  xlab("")


ggplot(filter(df), aes(x = time, group = scenario)) +
  geom_line(aes(y = LPI_boot, col = scenario)) +
  geom_line(aes(y = LPI_true, col = scenario), lty = 2) +
  facet_grid(Lag ~ direction)
             