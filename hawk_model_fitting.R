source("fitting_functions.R")

# Fit the c = 1 and c = 2 models to a sample ferruginous hawk migration (specifically, the migration displayed in Figure 1 of the paper). We also provide code for obtaining confidence intervals for the maximum likelihood estimates for each parameter.

require(tidyverse)

BUFFER = 7 # in days

hawk_data = read.csv("hawk_data_1hr.csv") # can also change filename to fit 12-hour and 24-hour timings
# Retain only the columns we need (x and y coordinates in UTM, plus previous coordinates for calculating step lengths and turn angles, plus the time)
hawk_data = as.matrix(hawk_data[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "dt")])

# Make the initial grid (with 14-day intervals) for possible t1 and t2 values to fit. Note that this looks different depending on c, the number of migrations
times_matrix_c1 = make_times_matrix(c(seq(min(hawk_data[, 7]), max(hawk_data[, 7]), by = 14), max(hawk_data[,7])), n_migrations = 1, min_diff = BUFFER)
times_matrix_c2 = make_times_matrix(c(seq(min(hawk_data[, 7]), max(hawk_data[, 7]), by = 14), max(hawk_data[,7])), n_migrations = 2, min_diff = BUFFER)

# Fit the model for c = 1 and c = 2. See the documentation within fitting_functions.R for more information on what each argument means
model_fit_c1 = estimate_times_nlevels(times_matrix_c1, c(14, 7, 3, 1), hawk_data, NLL_include = 5, bounds = c(min(hawk_data[, 7]), max(hawk_data[, 7])), min_diff = BUFFER)
model_fit_c2 = estimate_times_nlevels(times_matrix_c2, c(14, 7, 3, 1), hawk_data, NLL_include = 5, bounds = c(min(hawk_data[, 7]), max(hawk_data[, 7])), min_diff = BUFFER)

# The output contains parameter estimates for each point on the final subset of the grid. The true maximum likelihood estimate can be identified by finding the row with the lowest value in the "NLL" column.
MLE_c1 = model_fit_c1[which.min(model_fit_c1$NLL), ]
MLE_c2 = model_fit_c2[which.min(model_fit_c2$NLL), ]

# Obtain confidence intervals for the c = 1 data (takes a bit longer with c = 2 so I've opted to only show c = 1 here, but should be pretty easy to adjust the code if necessary).

N_BOOTSTRAP = 100 # how many simulations to use

best_par_indices = c(1, 2, 4, 5, 6, 7)
best_pars = as.numeric(MLE_c1[, best_par_indices]) # these are the indices of our parameters
names(best_pars) = colnames(MLE_c1)[best_par_indices]

# Bootstrapping. Note here that the "interval" argument includes more specific times to more precisely estimate the confidence bounds.
CIs_c1 = bootstrap_CI(hawk_data, best_pars, times_matrix_c1, n_iter = N_BOOTSTRAP, confidence_level = 0.95, intervals = c(14, 7, 3, 1, 1/2, 1/4, 1/12), NLL_include = 5)
