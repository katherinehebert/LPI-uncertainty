# check the propagated error confidence intervals

# load results
df <- readRDS("~/Documents/GitHub/LPI-sensitivity/outputs/all_results.RDS")

df$Process_error <- factor(df$Process_error, levels = c("0", "0.1", "0.2"))
df$interaction <- gsub("-0.2", "Strong Synchrony", df$interaction) 
df$interaction <- gsub("-0.1", "Weak Synchrony", df$interaction) 
df$interaction <- gsub("0.1", "Weak Asynchrony", df$interaction) 
df$interaction <- gsub("0.2", "Strong Asynchrony", df$interaction) 
df$interaction <- gsub("0", "No Synchrony", df$interaction) 

df$interaction <- factor(df$interaction, 
                         levels = c("Strong Asynchrony", "Weak Asynchrony", "No Synchrony", "Weak Synchrony", "Strong Synchrony"))



# set ggplot theme
theme_set(theme_linedraw())
format_lpiplots <- list(
  theme(legend.position = "none"),
  labs(x = "", y = "LPI"),
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 2)),
  ylim(c(0,1.8)))

# colour palette for interactions
colours <- c("Strong Asynchrony" = "#e66101", 
             "Weak Asynchrony" = "#fdb863", 
             "Weak Synchrony"= "#b2abd2", 
             "Strong Synchrony" = "#5e3c99")

# colour palette for process error plots
pal_lags <- rev(c("#22A884FF", "#2A788EFF", "#414487FF"))
pal_error <- c("#08306b", "#2171b5", "#6baed6")

# this was scenario 1, so filter for just that one
df_sc1 <- dplyr::filter(df, scenario %in% paste0("scenario1", LETTERS[1:9]))
(d <- ggplot(dplyr::filter(df_sc1, Process_error == "0"),
             aes(x = time, col = direction, group = scenario)) +
    geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se, 
                    fill = direction), alpha = .3, lwd = 0) +
    geom_line(aes(y = LPI_true), lty = 2, lwd = .2) +
    geom_line(aes(y = LPI_se)))
(d <- ggplot(dplyr::filter(df_sc1, Process_error == "0"),
             aes(x = time, col = direction, group = scenario)) +
    geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot, 
                    fill = direction), alpha = .3, lwd = 0) +
    geom_line(aes(y = LPI_true), lty = 2, lwd = .2) +
    geom_line(aes(y = LPI_boot)))


df_sc2 <- df[grepl("^scenario2", df$scenario), ] %>% dplyr::filter(Process_error == "0")
df_sc3 <- df[grepl("^scenario3", df$scenario), ] %>% dplyr::filter(Process_error == "0")
a <- ggplot(df_sc2,
            aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_boot)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_wrap(~interaction) +
  theme_bw() + format_lpiplots +
  labs(title = "Predation (negative covariation)")
b <- ggplot(df_sc3,
            aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_boot)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_wrap(~interaction) +
  theme_bw() + format_lpiplots +
  labs(title = "Competition (positive covariation)")
a + b

a <- ggplot(df_sc2,
            aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_se)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_wrap(~interaction) +
  theme_bw() + format_lpiplots +
  labs(title = "Predation (negative covariation)")
b <- ggplot(df_sc3,
            aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_se)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_wrap(~interaction) +
  theme_bw() + format_lpiplots +
  labs(title = "Competition (positive covariation)")
a + b

# make a lag df
df_lag <- dplyr::filter(df, 
              interaction %in% c("Weak Synchrony", "Strong Synchrony") & Lag %in% c("1", "2") & Process_error == "0")
ggplot(df_lag,
       aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_boot, ymax = cihi_boot, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_boot)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_grid(interaction~Lag, labeller = label_both) +
  theme_bw() + format_lpiplots
ggplot(df_lag,
       aes(x = time, group = scenario, col = direction)) +
  geom_ribbon(aes(ymin = cilo_se, ymax = cihi_se, 
                  fill = direction), alpha = .3, lwd = 0) +
  geom_line(aes(y = LPI_se)) +
  geom_line(aes(y = LPI_true), lty = 2) +
  facet_grid(interaction~Lag, labeller = label_both) +
  theme_bw() + format_lpiplots

ggplot(df) +
  geom_ribbon(aes(x = time, ymin = cilo_boot, ymax = cihi_boot, fill = scenario), alpha = .3) +
  geom_line(aes(x = time, y = LPI_boot, col = scenario)) +
  facet_wrap(direction ~ Process_error) +
  theme(legend.position = "none")
ggplot(df) +
  geom_ribbon(aes(x = time, ymin = cilo_boot, ymax = cihi_boot, fill = scenario), alpha = .3) +
  geom_line(aes(x = time, y = LPI_boot, col = scenario)) +
  facet_wrap(direction ~ interaction) +
  theme(legend.position = "none")
