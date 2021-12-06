# Script to gather all results into one table for plotting

# load packages and prep environment
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
theme_set(ggpubr::theme_pubr())

# loading outputs ----

# ## MY LPI ##
# 
# # import results of each scenario's LPI
# lpi <- lapply(paste0("outputs/", list.files(path = "outputs/", pattern = "_lpi.RDS")), readRDS)
# names(lpi) <- gsub("_lpi.RDS", "", list.files(path = "outputs/", pattern = "_lpi.RDS"))
# # and bind into one data frame
# lpi <- dplyr::bind_rows(lpi, .id = "scenario")

## RLPI ##

# import results from rlpi package
temp = list.files(path = "outputs/", pattern = "_rlpi.RDS")[-grep("true", list.files(path = "outputs/", pattern = "_rlpi.RDS"))]
rLpi <- lapply(paste0("outputs/", temp), readRDS)
names(rLpi) <- gsub("_rlpi.RDS", "", temp)
# and bind into one data frame
rLpi <- dplyr::bind_rows(rLpi, .id = "scenario")
rLpi$time <- rLpi$time - 1969

# import truth results from rlpi package
temp = list.files(path = "outputs/", pattern = "_rlpi.RDS")[grep("true", list.files(path = "outputs/", pattern = "_rlpi.RDS"))]
rLpi_true <- lapply(paste0("outputs/", temp), readRDS)
names(rLpi_true) <- gsub("_true_rlpi.RDS", "", temp)
# and bind into one data frame
rLpi_true <- dplyr::bind_rows(rLpi_true, .id = "scenario")
rLpi_true$time <- rLpi_true$time - 1969
colnames(rLpi_true)[2:4] <- paste0(colnames(rLpi_true)[2:4], "_true")

# join the rlpi results together
rLpi = right_join(rLpi, rLpi_true)

## PRECISION ####

temp = list.files(path = "outputs/", pattern = "_precision.RDS")
precision <- lapply(paste0("outputs/", temp), readRDS)
names(precision) <- gsub("_precision.RDS", "", temp)
# and bind into one data frame
precision <- dplyr::bind_rows(precision, .id = "scenario")
precision$time <- as.numeric(precision$time) - 1969
# join to results dataframe
rLpi = left_join(rLpi, precision)

## PARAMS ##

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
params$Lag[which(params$scenario %like% 'scenario4|scenario5')] <- "1"
params$Lag[which(params$scenario %like% 'scenario6|scenario7')] <- "2"

# join all tables together
df <- dplyr::left_join(rLpi, K_scenarios) %>% dplyr::left_join(params)

# accuracy ----

# calculate LPI accuracy 
df$accuracy <- df$LPI_final - df$LPI_final_true

# precision ----

# calculated in temporary.R

# format df columns for plotting
df$Lag <- factor(df$Lag, levels = c("0", "1", "2"))
df$direction <- factor(df$direction, levels = c("decline", "stable", "growth"))
colnames(df)[33] <- "N0"
colnames(df)[34] <- "lambda"
colnames(df)[35] <- "interaction"

# save to file
saveRDS(df, "outputs/all_results.RDS")
readr::write_csv(df, "outputs/all_results.csv")
