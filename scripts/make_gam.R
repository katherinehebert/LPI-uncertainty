# Function to build a generalized additive model (GAMs) per population

# load required packages
require(mgcv)
require(tidyverse)
require(errors)

# here we go!
make_gam <- function(simname){
  
  # import simulated populations
  pops <- readRDS(paste0("simulations/", simname, "_l.RDS")) %>%
    mutate_at(vars(time), as.integer)
  
  # 1. GAMs ----------------------------------------------------------------------
  
  # estimate the smoothing parameter to be half the time series length
  smoothParm = round(length(unique(pops$time))/2)
  
  # transform population sizes
  pops$N <- log10(pops$N + 1)
  # wrangle into wide format
  pops_w = pivot_wider(data = pops, id_cols = time, names_from = popID, values_from = N) %>%
    mutate_at(vars(time), as.integer) %>% as.matrix()
  
  # run GAM (modified code from CalcLPI function in rlpi package)
  m <- list()
  for(i in 2:ncol(pops_w)){
    N <- as.vector(pops_w[,i])
    time <- as.vector(pops_w[,"time"])
    m[[i-1]] <- gam(N ~ s(time, k = smoothParm), 
                    family = gaussian(), fx = TRUE, method = "REML")
  }
  names(m) <- unique(pops$popID)
  saveRDS(m, paste0("models/", simname, "_gam.RDS"))
  
  # predict over time period
  pred_ls = lapply(m, predict.gam, type = "response", se.fit = TRUE)
  
  # wrangle into one long format dataframe
  pred = pred_ls %>%
    lapply(bind_cols) %>%
    lapply(dplyr::mutate, time = pops_w[,"time"]) %>%
    bind_rows(.id = "popID")
  # join to observation dataframe
  pred = full_join(pops, pred, by = c("popID", "time"))
  pred = dplyr::rename(pred,
                "N_pred" = "fit",
                "N_se" = "se.fit")
  
  # plot the predicted trend over the original
  ggplot(pred, aes(x = time, group = popID)) +
    #geom_line(aes(y = N, col = pop), lty = 2, lwd = .5) +
    geom_ribbon(aes(ymin = N_pred - N_se,
                    ymax = N_pred + N_se, fill = pop), alpha = .1) +
    geom_line(aes(y = N_pred, col = pop)) +
    labs(y = "N(t) (log10)", col = "Populations", fill = "Populations",
         caption = "Ribbon shows standard error from GAM predictions.") + 
    facet_wrap(~set) +
    theme(legend.position = "none")
  ggsave(filename = paste0(simname, "_gampred.png"), path = "figures/", plot = last_plot(),
         width = 7, height = 5, units = "in")
  
  # calculate dt from GAM predictions
  dt <- lapply(pred_ls, get_dt, time = pops_w[,"time"]) %>% bind_rows(.id = "popID")
  # join to observations dataframe
  dt_df = full_join(pred, dt, by = c("popID", "time"))
  
  # plot dts
  ggplot(data = dt_df, aes(x = time, group = popID)) +
    geom_ribbon(aes(ymin = dt-se, ymax = dt+se, fill = pop), alpha = .1) +
    geom_line(aes(y=dt, col = pop)) + 
    facet_wrap(~set) + 
    labs(y = "Growth rate (log10)", 
         col = "Populations", fill = "Populations",
         caption = "Ribbon shows propagated standard error from GAM predictions.") +
    theme(legend.position = "none")
  
  saveRDS(dt_df, paste0("outputs/", simname, "_results.RDS"))
  ggsave(filename = paste0(simname, "_dt.png"), path = "figures/", plot = last_plot(),
         width = 7, height = 5, units = "in")
}


