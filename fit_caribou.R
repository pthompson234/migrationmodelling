require(Rcpp)

sourceCpp("cpp_functions.cpp")
source("fitting_functions.R")

N_MIGRATIONS = 1
BUFFER = 14

data_in = read.csv("sample_caribou.csv")[, -1]

data_for_migration = as.matrix(data_in[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "JUL")])
max_time = min(max(data_for_migration[,7]), 180) # force it to end at the end of June at the latest because we're looking only for the spring migration
data_for_migration = data_for_migration[data_for_migration[,7] <= max_time, ]

# For traditional fits
times_matrix_free = make_times_matrix(starts = c(seq(min(data_for_migration[, 7]), max(data_for_migration[, 7]), by = 14), max(data_for_migration[, 7])), n_migrations = N_MIGRATIONS, min_diff = BUFFER)
model_fit_free = estimate_times_nlevels(times_matrix_free, c(14, 7, 3, 1), mod_type = "r1k1", data_for_migration, NLL_include = 5, cpp = TRUE, bounds = c(min(data_for_migration[, 7]), max(data_for_migration[, 7])), min_diff = BUFFER)

# For fits with t2 fixed (we know the calving date / end of migration already)
this_calv_date = 155

times_matrix_fixed = make_times_matrix(starts = c(seq(min(data_for_migration[, 7]), max(data_for_migration[, 7]), by = 1), max(data_for_migration[,7])), ends = this_calv_date, n_migrations = N_MIGRATIONS, min_diff = BUFFER)
model_fit_fixed = estimate_times_cpp(times_matrix_fixed, data_for_migration)

model_fit_free = model_fit_free[which.min(model_fit_free[, 3]), , drop = FALSE]
model_fit_fixed = model_fit_fixed[which.min(model_fit_fixed[, 3]), , drop = FALSE]

colnames(model_fit_free) = c("t1", "t2", "NLL", "rho1", "rho2", "kappa1", "kappa2")
colnames(model_fit_fixed) = names(model_fit_free)