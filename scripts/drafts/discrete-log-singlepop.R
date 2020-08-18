# Discrete logistic population growth simulation

# parameters for simulation ----------------------------------------------------

# max growth rate
r = 1
# current population size
N0 = 10
# carrying capacity
K = 500
# number of generations to simulate
gen = 10


# run simulation ---------------------------------------------------------------

# initialize vector to store results (population sizes)
N <- N0
# initialize variable to hold current population size
N_t <- N0

# calculate population sizes
for(i in 1:gen){
  N_t <- N_t + N_t*r*(1-N_t/K)
  if(N_t < 0) {N_t <- 0} # assign 0 to negative population sizes
  N <- c(N, N_t) # append resulting population size to results vector
}

# plot results -----------------------------------------------------------------

time <- 1:length(N)
plot(N ~ time, 
     type = "o", pch = 16, las = 1,
     xlab = "Time", ylab = "Population size (N)")
