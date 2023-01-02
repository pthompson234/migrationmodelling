source("fitting_functions.R")

# Simulate a migration (as was done in the paper) and fit the migration model (with c = 1) to the simulated data.

require(adehabitatLT)
require(changepoint)
require(bcpa)
require(tidyverse)

set.seed(123) # Note this will produce different results from the actual paper but will provide consistency for those using the Github repo.

# True values for all model parameters
t1 = 70
t2 = 100

N_DAYS = 200
N_STEPS_DESIRED = 150

RHO_0 = 5
RHO_1 = 40
KAPPA_0 = 0
KAPPA_1 = 0.5

BUFFER = 7
FPT_RADII = seq(5, 100, 5) # set of first passage time radii to test

# Simulate random data
simulated_data = sim_migratory_movement(N_DAYS, t1, t2, percent_removed = (1 - N_STEPS_DESIRED / N_DAYS) / 2, sl_migration = RHO_0 + RHO_1, sl_resident = RHO_0, kappa_migration = KAPPA_1, kappa_resident = KAPPA_0, return_code2_only = TRUE)
simulated_data = cbind(simulated_data, Code = 2)

# Fit our model
times_matrix = make_times_matrix(c(seq(min(simulated_data[, 7]), max(simulated_data[, 7]), by = 14), max(simulated_data[,7])), n_migrations = 1, min_diff = BUFFER)
model_fit_thompson = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), simulated_data, NLL_include = 5, bounds = c(min(simulated_data[, 7]), max(simulated_data[, 7])), min_diff = BUFFER)
model_fit_thompson = model_fit_thompson[which.min(model_fit_thompson$NLL), ]

# Fit the NSD model from Bunnefeld et al (2011)
model_fit_bunnefeld = fit_NSD_migration(simulated_data[, c(1,2,7)], units = "km")

# Fit the FPT model from Le Corre et al. (2014)
model_fit_lecorre = fit_fpt_penalized_constant(simulated_data[, c(1,2,7)], FPT_RADII)

# Fit the PELT model using daily distances from Madon & Hingrat (2014)
model_fit_madon = fit_dailydist_pelt(simulated_data)

# Feel free to compare the model fits, particularly the estimates of t1 and t2 from each model.