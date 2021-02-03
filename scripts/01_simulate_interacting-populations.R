# Script to simulate population size over time for multiple paired populations 
# that are interacting together

# Interaction coefficient is constant through time
# All populations have the same "true" growth rate
# Observation error is added to the growth rate (varies randomly through time)
# Process error is added to the current population size (varies randomly thorugh time)

# load packages
require(tidyverse)
# set ggplot theme
theme_set(theme_linedraw() + theme(panel.grid = element_blank()))
# set seed for randomization
set.seed(2)

# parameters for simulation ----------------------------------------------------

# number of population pairs to generate
n_pairs = 10

# number of generations
gen = 10 

# initial population sizes
N0i = 200
N0j = 200

# growth rate
lambda_i = 0.5 # true growth rate
lambda_j = 0.5 # true growth rate
# matrix form
r = matrix(lambda_i, nrow = n_pairs, ncol = gen) 

# interaction coefficients (competition)
alpha_ij = -0.1
alpha_ji = -0.2

# carrying capacity
Ki = 190
Kj = 190

## ERROR ##

process = runif(n_pairs*gen, -N0i/20, N0i/20) # process error (environmental noise)
observation = runif(n_pairs*gen, 0, 0) # observation error

# add process error to growth rates
r_error = r + process
# observation error
obs_error = matrix(observation, nrow = n_pairs, ncol = gen)

# run simulation --------------------------------------------------------------- 

# initialize matrix to store results (population sizes)
Ni <- as.matrix(rep(N0i, n_pairs))
Nj <- as.matrix(rep(N0j, n_pairs))

# Nt+1_i = Nt_i + rNt_i * ((1 - Nt_i/K_i) + alpha_ji*Nt_j/K_j
# calculate population sizes
for(t in 1:gen-1){
  
  # population i
  temp_i = Ni[t]*(1 + r_error[,t]*(1 - (Ni[t] + alpha_ij*Nj[t])/Ki)) + obs_error[,t]
  #temp_i[which(temp_i < 0)] <- 0 # assign 0 to negative population sizes
  Ni <- cbind(Ni, temp_i) # append resulting population size to results vector
  
  # population j
  temp_j = Nj[t]*(1 + r_error[,t]*(1 - (Nj[t] + alpha_ji*Ni[t])/Kj)) + obs_error[,t]
  #temp_j[which(temp_j < 0)] <- 0 # assign 0 to negative population sizes
  Nj <- cbind(Nj, temp_j) # append resulting population size to results vector
}

# plot results -----------------------------------------------------------------

# create vector of time values for plotting
time <- 1:gen

# function to wrangle the results into long format 
pops_long <- function(pops_df, n = n_pairs, g = gen, set_id) {
  pops_df = as.data.frame(pops_df)
  colnames(pops_df) = time
  pops_df = mutate(.data = pops_df, "popID" = paste(set_id, sprintf("pop%s", 1:n), sep = "-")) %>%
    pivot_longer(cols = 1:all_of(g), names_to = "time", values_to = "N") %>%
    separate(popID, into = c("set", "pop"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(time), as.integer)
}

# bind together
N = rbind(pops_long(Ni, set_id = "i"), pops_long(Nj, set_id = "j"))

# plot
ggplot(N) +
  geom_line(aes(x = time, y = N, group = popID, col = popID)) + 
  #ylim(c(min(N$N - 10), max(N$N+10))) +
  #geom_hline(yintercept = Ki, lty = 2) +
  #geom_hline(yintercept = Kj, lty = 2) +
  facet_wrap(~ set) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(N$N)+10))

# save outputs -----------------------------------------------------------------
saveRDS(N, "simulations/paired_antagonistic_l.RDS")
ggsave(filename = "paired_antagonistic_N.png", path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")

# calculate covariation --------------------------------------------------------

N_w = cbind(t(Ni), t(Nj))

# plot covariation
png("figures/paired_antagonistic_cov.png", width = 500, height = 500)
cov(N_w) %>% heatmap(Colv = NA, Rowv = NA, 
                        col = (colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)),
                        main = "Covariation")
dev.off()

# plot correlation
png("figures/paired_antagonistic_cor.png", width = 500, height = 500)
cor(N_w) %>% heatmap(Colv = NA, Rowv = NA, 
                        col = (colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)), 
                        main = "Correlation")
dev.off()