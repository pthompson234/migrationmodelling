source("fitting_functions.R")

# Fit the migration model to the caribou data depicted in Figure 5 of the paper. We compare our results to the NSD model and also use a second times_matrix where t2 is fixed, based on the calving date which we identify using a simple broken-stick regression model (DeMars et al., 2013).

caribou_data = read.csv("caribou_data.csv")
caribou_data = as.matrix(caribou_data[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "JUL")])

max_time = min(max(caribou_data[,7]), 180) # Only estimate t1 and t2 between Jan 1 and Jun 30 as we are focused on spring migration

# Fit the "traditional" c = 1 version of the model by estimating both t1 and t2.
times_matrix_orig = make_times_matrix(starts = c(seq(min(caribou_data[, 7]), max_time, by = 14), max_time), n_migrations = 1)
model_fit_orig = estimate_times_nlevels(times_matrix_orig, c(14, 7, 3, 1), caribou_data, NLL_include = 5, bounds = c(min(caribou_data[, 7]), max(caribou_data[, 7])))

# Fit the NSD model to the data
model_fit_NSD = fit_NSD_migration(caribou_data[caribou_data[,7] < max_time, c(1,2,7)], units = "m")

# Obtain the calving date
test_times = make_times_matrix(starts = 145:180, ends = 150:240, min_diff = 7)
calving_date_estimates = fit_calving(caribou_data, test_times)
calving_date = calving_date_estimates$t1[which.min(calving_date_estimates$NLL)][[1]]

# Fit the c = 1 model again but this time, with a new times_matrix fixing t2 at the calving date
times_matrix_t2fixed = make_times_matrix(starts = c(seq(min(caribou_data[, 7]), max(caribou_data[, 7]), by = 1), max(caribou_data[,7])), ends = calving_date, n_migrations = 1)
model_fit_t2fixed = estimate_times(times_matrix_t2fixed, caribou_data[caribou_data[,7] <= calving_date, ])