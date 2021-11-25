# script to check that the calculated LPI corresponds to the same as would be
# obtained from the rlpi package

# i.e. calculate the LPI with the rlpi functions to ensure it's correct

# my CI are bigger than the rlpi's
# mine is slighly higher when it's decline, and slightly lower when it's growth

# so there must be rounding issues, but it is still quite similar

# clear workspace
rm(list=ls())

# Load library
library(rlpi)
library(dplyr)
library(tidyr)
library(ggplot2)

# source customised LPI function
source('~/Documents/GitHub/LPI-sensitivity/scripts/LPIMain.R')

theme_set(theme_linedraw())

lpi_checker <- function(scenario_name){
  
  # import simulated population
  sim <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/simulations/", scenario_name, "_l.RDS"))
  
  # change time column to fake years, for the rlpi function to work
  sim$time <- paste0("X", sim$time + 1969)
  
  # convert to wide format
  simw <- pivot_wider(sim, names_from = time, values_from = N)
  
  # create ID and Binomial columns
  simw$ID <- 1:20
  simw$Binomial <- LETTERS[1:20]
  
  # Constructing infiles from a populations table
  
  # Select populations in the dataset by setting all to TRUE
  index_vector = rep(TRUE, nrow(simw))
  
  # Create infile for Canada populations
  checker_infile_name <- create_infile(simw, 
                                       index_vector = index_vector, 
                                       name = scenario_name, 
                                       start_col_name = "X1970", 
                                       end_col_name = "X1980"
                                       )
  
  
  # index with 1000 bootstraps (1137 pops) without weightings 
  lpi <- LPIMain_custom(checker_infile_name, 
                 REF_YEAR = 1970, PLOT_MAX = 1980, 
                 BOOT_STRAP_SIZE = 1000, 
                 use_weightings=0, 
                 VERBOSE=FALSE, 
                 save_plots = 0, 
                 plot_lpi = 0,
                 basedir = "outputs/rlpi",
                 scenario_name = scenario_name, 
                 force_recalculation = TRUE)
  # Remove NAs (trailing years with no data)
  lpi <- lpi[complete.cases(lpi), ]
  lpi$time <- rownames(lpi) %>% as.numeric()
  saveRDS(lpi, paste0("outputs/", scenario_name, "_rlpi.RDS"))

  # plot this LPI
  ggplot(data = lpi, aes(x = time)) +
    geom_ribbon(data = lpi, aes(ymin = CI_low, 
                                ymax = CI_high), alpha = .2, fill = "navyblue") +
    geom_line(data = lpi, aes(y = LPI_final), col = "navyblue") +
    labs(x = "", y = "LPI") +
    theme(legend.position = "none")
  ggsave(paste0("figures/", scenario_name, "_rlpi.png"))
  
  
  ## check the difference between LPI and my calculated version...
  
  sim_lpi <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/", scenario_name, "_lpi.RDS"))
  sim_lpi$time <- sim_lpi$time + 1969
  
  # Plot the resulting indices to compare them, and save it
  
  # make palette with names
  pal = c("rLPI" = "red", "myLPI" = "blue")
  ggplot() +
    # rlpi LPI
    geom_ribbon(data = lpi, aes(x = time, 
                                ymin = CI_low, 
                                ymax = CI_high,
                                fill = "rLPI"), alpha = .2) +
    geom_line(data = lpi, aes(x = time, 
                              y = LPI_final, 
                              col = "rLPI")) + 
    # my calculation
    geom_ribbon(data = sim_lpi, aes(x = time, 
                                    ymin = cilo_boot, 
                                    ymax = cihi_boot, 
                                    fill = "myLPI"), 
                alpha = .2) +
    geom_line(data = sim_lpi, aes(x = time, 
                                  y = LPI_boot, 
                                  col = "myLPI")) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    labs(fill = "Source", col = "Source")
  ggsave(paste0("figures/", scenario_name, "_lpicheck.png"))
  
  
  # make a table with the results
  checker_df <- data.frame(
    "scenario" = scenario_name,
    "LPI_diff" = sim_lpi$LPI_boot - lpi$LPI_final,
    "cilo_diff" = sim_lpi$cilo_boot - lpi$CI_low,
    "cihi_diff" = sim_lpi$cihi_boot - lpi$CI_high,
    "time" = sim_lpi$time - 1969
  )
  saveRDS(checker_df, paste0("outputs/", scenario_name, "_lpicheck.RDS"))
  
}


# get all scenario names
sim_names <- c(
  paste0("scenario1", LETTERS[1:9]),
  paste0("scenario2", LETTERS[1:18]),
  paste0("scenario3", LETTERS[1:18]),
  paste0("scenario4", LETTERS[1:18]),
  paste0("scenario5", LETTERS[1:18]),
  paste0("scenario6", LETTERS[1:18]),
  paste0("scenario7", LETTERS[1:18])
  )
# apply to all simulations
lapply(sim_names, lpi_checker)              

               
### calculate true LPI with rlpi

lpi_truechecker <- function(scenario_name){
  
  # import simulated population
  sim <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/", scenario_name, "_true.RDS"))
  
  # change time column to fake years, for the rlpi function to work
  sim$time <- paste0("X", sim$time + 1969)
  
  # convert to wide format
  sim <- subset(sim, select = -dt)
  simw <- pivot_wider(sim, names_from = time, values_from = N)
  
  # create ID and Binomial columns
  simw$ID <- 1:20
  simw$Binomial <- LETTERS[1:20]
  
  # Constructing infiles from a populations table
  
  # Select populations in the dataset by setting all to TRUE
  index_vector = rep(TRUE, nrow(simw))
  
  # Create infile for Canada populations
  checker_infile_name <- create_infile(simw, 
                                       index_vector = index_vector, 
                                       name = paste0(scenario_name, "_true"), 
                                       start_col_name = "X1970", 
                                       end_col_name = "X1980"
                                       )
  
  
  # index with 1000 bootstraps (1137 pops) without weightings 
  lpi <- LPIMain(checker_infile_name, 
                 REF_YEAR = 1970, PLOT_MAX = 1980, 
                 BOOT_STRAP_SIZE = 1000, 
                 use_weightings=0, 
                 VERBOSE=FALSE, save_plots = 0, plot_lpi = 0,
                 basedir = "outputs/rlpi", force_recalculation = TRUE)
  # Remove NAs (trailing years with no data)
  lpi <- lpi[complete.cases(lpi), ]
  lpi$time <- rownames(lpi) %>% as.numeric()
  saveRDS(lpi, paste0("outputs/", scenario_name, "_true_rlpi.RDS"))
  
  # plot this LPI
  ggplot(data = lpi, aes(x = time)) +
    geom_ribbon(data = lpi, aes(ymin = CI_low, 
                                ymax = CI_high), alpha = .2, fill = "navyblue") +
    geom_line(data = lpi, aes(y = LPI_final), col = "navyblue") +
    labs(x = "", y = "LPI") +
    theme(legend.position = "none")
  ggsave(paste0("figures/", scenario_name, "_true_rlpi.png"))
  
  
  ## check the difference between LPI and my calculated version...
  
  sim_lpi <- readRDS(paste0("~/Documents/GitHub/LPI-sensitivity/outputs/", scenario_name, "_lpi.RDS"))
  sim_lpi$time <- sim_lpi$time + 1969
  
  # Plot the resulting indices to compare them, and save it
  
  # make palette with names
  pal = c("rLPI" = "red", "myLPI" = "blue")
  ggplot() +
    # rlpi LPI
    geom_line(data = lpi, aes(x = time, 
                              y = LPI_final, 
                              col = "rLPI")) + 
    # my calculation
    geom_line(data = sim_lpi, aes(x = time, 
                                  y = LPI_boot, 
                                  col = "myLPI")) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    labs(fill = "Source", col = "Source")
  ggsave(paste0("figures/", scenario_name, "_lpicheck_true.png"))
  
  
  # make a table with the results
  checker_df <- data.frame(
    "scenario" = scenario_name,
    "LPI_diff" = sim_lpi$LPI_boot - lpi$LPI_final,
    "time" = sim_lpi$time - 1969
  )
  saveRDS(checker_df, paste0("outputs/", scenario_name, "_lpicheck_true.RDS"))
  
}

# get all scenario names
sim_names <- c(
  paste0("scenario1", LETTERS[1:9]),
  paste0("scenario2", LETTERS[1:18]),
  paste0("scenario3", LETTERS[1:18])#,
  paste0("scenario4", LETTERS[1:18]),
  paste0("scenario5", LETTERS[1:18]),
  paste0("scenario6", LETTERS[1:18]),
  paste0("scenario7", LETTERS[1:18])
)
lapply(sim_names, lpi_truechecker)
