# Script to gather all results into one table for plotting

# load packages and prep environment
library(ggplot2)
library(dplyr)
library(data.table)
theme_set(ggpubr::theme_pubr())

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
    x <- x[-c(1,2),]
  })
names(params) <- gsub("_params.RDS", "", list.files(path = "simulations/", pattern = "_params.RDS"))
# and bind into one data frame
params <- dplyr::bind_rows(params, .id = "scenario")

# make table for categories from carrying capacity scenarios
K_scenarios <- data.frame("scenario" = params$scenario, "direction" = NA)
K_scenarios$direction[which(K_scenarios$scenario %like% 'A|D|G|J|M|P')] <- "decline"
K_scenarios$direction[which(K_scenarios$scenario %like% 'B|E|H|K|N|Q')] <- "stable"
K_scenarios$direction[which(K_scenarios$scenario %like% 'C|F|I|L|O|R')] <- "growth"

# associate lag with the scenarios
# params$Lag[which(params$scenario %like% 'scenario4|scenario5')] <- "1"
# params$Lag[which(params$scenario %like% 'scenario6|scenario7')] <- "2"

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
df$interval_width <- ((df$cihi_boot-df$cilo_boot)/df$LPI_boot)*100
df$interval_width_se <- ((df$cihi_se-df$cilo_se)/df$LPI_se)*100

# format columns for plotting
df$Lag <- factor(df$Lag, levels = c("0", "1", "2"))
df$direction <- factor(df$direction, levels = c("decline", "stable", "growth"))
colnames(df)[11] <- "N0"
colnames(df)[12] <- "lambda"
colnames(df)[13] <- "interaction"

### TEMPORARY until simulation is fixed for 3M and 5M onwards
df <- df[-which(df$scenario %in% c(paste0("scenario3", LETTERS[13:18]), paste0("scenario5", LETTERS[13:18]))),]

# save to file
saveRDS(df, "outputs/all_results.RDS")
readr::write_csv(df, "outputs/all_results.csv")

##### ALSO: compare interval width when error is propagated vs. when it's bootstrapped!!!! missing a bunch of error aren't we?