require(Rcpp)
require(ctmm)
require(tidyverse)

source("fitting_functions.R")
sourceCpp("cpp_functions.cpp")

N_MIGRATIONS = 2

data_in = read.csv("sample_bear.csv")[, -1]
data_for_migration = as.matrix(data_in[, c("easting", "northing", "xprev", "yprev", "xprevprev", "yprevprev", "jday")])

times_matrix = make_times_matrix(c(seq(min(data_for_migration[, 7]), max(data_for_migration[, 7]), by = 14), max(data_for_migration[,7])), n_migrations = N_MIGRATIONS)

model_fit = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), data_for_migration, mod_type = "r1k1", NLL_include = 5, cpp = TRUE, bounds = c(min(data_for_migration[, 7]), max(data_for_migration[, 7])))
model_fit = model_fit[which.min(model_fit[,5]), , drop = FALSE]
colnames(model_fit) = c("t1", "t2", "t3", "t4", "NLL", "r0", "r1", "k0", "k1")