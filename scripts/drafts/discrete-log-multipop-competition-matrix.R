# Multispecies discrete logistic population growth simulation with competition
# for multiple pairs of populations (matrix!)

# load packages
require(tidyverse)
# set ggplot theme
theme_set(theme_linedraw() + theme(panel.grid = element_blank()))
# set seed for randomizations
set.seed(2)

# parameters for simulation ----------------------------------------------------

n_pairs = 10 # number of population pairs to generate
gen = 10 # number of generations
N0 = 100
process = runif(n_pairs*gen, -N0/10, N0/10) # process error (environmental noise)
observation = runif(n_pairs*gen, -0.1, 0.1) # observation error

# current population sizes (t=0)
N0_A = rep(N0, n_pairs)
N0_B = rep(N0, n_pairs)

# growth rate
r = matrix(1, nrow = n_pairs, ncol = gen) 
# add observation error to growth rates
r_error = r + observation

# interaction coefficients (competition)
alpha_AB = -0.0002
alpha_BA = 0.0002

# process error (environmental noise)
proc_error = matrix(process, nrow = n_pairs, ncol = gen)

# run simulation --------------------------------------------------------------- 

# initialize matrix to store results (population sizes)
popsA <- as.matrix(N0_A)
popsB <- as.matrix(N0_B)

# calculate population sizes
# current population size * (growth rate + observation error) + environmental noise
for(t in 2:gen){
  # set A
  temp_A <- popsA[,t-1]*(r_error[,t] + alpha_BA*popsB[,t-1]) + proc_error[,t]
  temp_A[which(temp_A < 0)] <- 0 # assign 0 to negative population sizes
  popsA <- cbind(popsA, temp_A) # append resulting population size to results vector
  # set B
  temp_B <- popsB[,t-1]*(r_error[,t] + alpha_AB*popsA[,t-1]) + proc_error[,t]
  temp_B[which(temp_B < 0)] <- 0 # assign 0 to negative population sizes
  popsB <- cbind(popsB, temp_B) # append resulting population size to results vector
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
    separate(popID, into = c("set", "pop"), sep = "-", remove = FALSE)
}

# bind together
pops = rbind(pops_long(popsA, set_id = "A"), pops_long(popsB, set_id = "B"))

# plot
ggplot(pops) +
  geom_line(aes(x = time, y = N, group = popID, col = pop)) + 
  ylim(c(0, max(pops$N+10))) +
  facet_wrap(~set) + theme(legend.position = "none")

# save outputs -----------------------------------------------------------------
saveRDS(pops_long, "outputs/paired_antagonistic_l.RDS")
ggsave(filename = "paired_antagonistic_N.png", path = "figures/", plot = last_plot(),
       width = 7, height = 5, units = "in")
