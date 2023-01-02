source("fitting_functions.R")

# Fit the migration model to a movement track from a brown bear in the Mackenzie Delta. We use c = 2 and fit the model to the entire year of movement.

bear_data = read.csv("bear_data.csv")
bear_data = as.matrix(bear_data[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "jday")])

# Fit the c = 2 version of the model.
times_matrix = make_times_matrix(c(seq(min(bear_data[, 7]), max(bear_data[, 7]), by = 14), max(bear_data[,7])), n_migrations = 2)

model_fit = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), bear_data, NLL_include = 5, bounds = c(min(bear_data[, 7]), max(bear_data[, 7])))
