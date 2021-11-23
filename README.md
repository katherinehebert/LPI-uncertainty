# LPI-sensitivity

 Sensitivity tests on the Living Planet Index

Steps:

1. run_sims.R simulates population time series, calculates the LPI of these populations, and simulates a "true" set of populations without observation error.
2. check_lpi.R calculates the LPI of the simulated populations with error and the true populations using the rlpi package. I have found some differences between my LPI and this rlpi version, so I will be basing myself on the rlpi results.
3. gather_results.R gathers the results of these scripts into a dataframe of results about the accuracy and precision of the different LPIs.
4. get_chaindt.R calculates the log-ratio growth rates of the simulated (error) populations for some comparisons.
5. docs/accuracy_and_precision.Rmd then plots these results in various ways in their final form, which can then be used in the manuscript (in Ch1_manuscript repo).

To do:

temporary.R has code near the bottom to calculate the precision of the LPI, and it needs to be cleaned up into its own script and integrated into gather_results.R. 
