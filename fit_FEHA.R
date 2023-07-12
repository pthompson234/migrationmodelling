require(Rcpp)
require(tidyverse)

sourceCpp("cpp_functions.cpp")
debugSource("fitting_functions.R")

BUFFER = 7
N_MIGRATIONS = 1 # Feel free to change this to 2 and see what happens

data_in = read_csv("sample_FEHA.csv")[, -1]
data_for_migration = as.matrix(data_in[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "dt")])

# This is for n_levels
times_matrix = make_times_matrix(c(seq(min(data_for_migration[, 7]), max(data_for_migration[, 7]), by = 14), max(data_for_migration[,7])), n_migrations = N_MIGRATIONS, min_diff = BUFFER)

# Model fitting the new way (with n levels)
model_fit = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), data_for_migration, mod_type = "r1k1", NLL_include = 5, bounds = c(min(data_for_migration[, 7]), max(data_for_migration[, 7])), min_diff = BUFFER, cpp = TRUE)
this_migration = model_fit[which.min(model_fit[,1+2*N_MIGRATIONS]), , drop = FALSE]
best_pars = as.numeric(this_migration[, -(1 + ncol(times_matrix))]) # these are the indices of our parameters (can check with the all_fits file)
names(best_pars) = c(paste0("t", 1:(2*N_MIGRATIONS)), "r0", "r1", "k0", "k1")

# If N_MIGRATIONS is 2 this takes a few minutes.
CIs = bootstrap_CI(data_in = data_for_migration, optim_pars = best_pars, times_matrix = times_matrix, n_iter = 100, intervals = c(14, 7, 3, 1), NLL_include = 5, cpp = TRUE)
