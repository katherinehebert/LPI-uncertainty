# Testing the Living Planet Index's ability to capture complex biodiversity change

Code for simulations, calculations, and figures presented in the manuscript "Testing the Living Planet Index's ability to capture complex biodiversity change" by Hébert, Katherine & Gravel, Dominique.

To produce the results in the manuscript, these scripts can be run:

`01_run_sims.R`: Simulates population time series and simulates a "true" set of populations without observation error.

`02_calculate_lpi.R`: Calculates the LPI of the simulated populations with error and the true populations using the rlpi package. 

`03_get_chaindt.R`: Calculates the log-ratio growth rates of the simulated populations for some comparisons with the smoothed growth rates in the LPI calculation.

`04_propagate_uncertainty.R`: Calculates the propagated uncertainty from noise in the raw data, through to the final LPI trend, following equations developed in the manuscript.

`05_get_precision.R`: Calculates some measures of precision about the LPI trend.

`06_gather_results.R`: Gathers the results of these scripts into a dataframe of results about the accuracy and uncertainty of the LPIs, in preparation for plotting.

`07_make_figures.R`: Generates the figures presented in the main manuscript, and the supplementary matrials.
