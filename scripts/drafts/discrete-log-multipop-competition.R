# Multispecies discrete logistic population growth simulation with competition

# parameters for simulation ----------------------------------------------------

# max growth rate
r_A = .1
r_B = .1

# current population size
N0_A = 1000
N0_B = 1000

# carrying capacity
K_A = 500
K_B = 500

# interaction coefficients (competition)
alpha_A = .1
alpha_B = .9

# number of generations to simulate
gen = 10


# run simulation ---------------------------------------------------------------

# initialize vector to store results (population sizes)
popA <- N0_A
popB <- N0_B

# initialize variable to hold current population size
N_A <- N0_A
N_B <- N0_B

# calculate population sizes
for(i in 1:gen){
  
  # population A
  N_A <- N_A + N_A*r_A*(1-N_A/K_A - alpha_A*N_B/K_A)
  if(N_A < 0) {N_A <- 0} # assign 0 to negative population sizes
  popA <- c(popA, N_A) # append resulting population size to results vector
  
  # population B
  N_B <- N_B + N_B*r_B*(1-N_B/K_B - alpha_B*N_A/K_B)
  if(N_B < 0) {N_B <- 0} # assign 0 to negative population sizes
  popB <- c(popB, N_B) # append resulting population size to results vector
  
}



# plot results -----------------------------------------------------------------

#find maximum value for y axis
ymax <- max(c(max(popA), max(popB)))
YAxisLength <- max(pretty(ymax*1.05))

# create vector of time values for plotting
time <- 1:length(popA)

# plot
plot(popA ~ time, 
     type = "o", col = "red", ylim = c(0,YAxisLength),
     xlab = "Time", ylab = "Population size (N")
points(popB ~ time, col = "blue", type = "o")
abline(h = K_A, lty = 2, col = "red")
abline(h = K_B, lty = 2, col = "blue")