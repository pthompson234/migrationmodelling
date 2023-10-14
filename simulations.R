library(adehabitatLT)
library(changepoint)
library(bcpa)
library(tidyverse)
library(mcp)
library(marcher)
library(rjags)
library(Rcpp)

source("fitting_functions.R")
sourceCpp("cpp_functions.cpp")
# functions from waddle package
waddle_dir = "" # replace with where the waddle library is located in your device 
lapply(list.files(waddle_dir), function(f) {
  source(paste0(waddle_dir, "/R/", f))
})

set.seed(111) # for consistency

N_SIM = 50 # How many random migrations to simulate? In reality we want 50; this is coverage for the occasional failed convergence from the NSD model.

N_STEPS = 300 # number of timesteps
T1 = 100 # beginning of migration(s)
T2 = 200 # end of migration(s)
PERCENT_REMOVED = 1/12 # gives us about 250 steps
ERROR = 0 # gps error; let's start by assuming it's negligible and move from there

# Which model are we simulating from? Pick the simulation function here and simulate from one of the 3 processes from Gurarie et al 2016
sim_model = "speed"

all_results = data.frame()

# TO DO: 1) Is FPT working the way it's supposed to? Test one of the bias fits (in Gurarie et al 2016 it did great with those; maybe they do it differently than I did or maybe mine has a bug?)

for (i in 1:N_SIM) {
  
  # simulate the path; depends on what model we choose. Speciic parameters taken roughly from Gurarie et al 2016
  if (sim_model == "speed") {
    data_for_migration = sim_migration_cvm(n_timesteps = N_STEPS, migration_starts = T1, migration_ends = T2, percent_removed = PERCENT_REMOVED, error = ERROR, nu_resident = 1, nu_migration = 5, tau_resident = 2, tau_migration = 2, return_code2_only = FALSE)
  } else if (sim_model == "timescale") {
    data_for_migration = sim_migration_cvm(n_timesteps = N_STEPS, migration_starts = T1, migration_ends = T2, percent_removed = PERCENT_REMOVED, error = ERROR, nu_resident = 1, nu_migration = 1, tau_resident = 2, tau_migration = 20, return_code2_only = FALSE)
  } else { # if (sim_model == "bias")
    data_for_migration = sim_migration_bcrw(n_timesteps = N_STEPS, migration_starts = T1, migration_ends = T2, migration_distances = 50, percent_removed = PERCENT_REMOVED, error = ERROR, sl_migration = 1, att_migration = 0.9, att_resident = 0.5, return_code2_only = FALSE)
  }
  
  t0 = Sys.time()
  
  # Our model
  data_code2only = data_for_migration %>% dplyr::filter(Code == 2) %>% dplyr::select(-Code) %>% as.matrix
  times_matrix = make_times_matrix(c(seq(min(data_code2only[, 7]), max(data_code2only[, 7]), by = 14), max(data_code2only[,7])), n_migrations = 1, min_diff = 7)
  model_fit_mmcp = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), data_code2only, 
                                          mod_type = "r1k1", NLL_include = 5, cpp = TRUE, min_diff = 7,
                                          bounds = c(min(data_code2only[, 7]), max(data_code2only[, 7])))
  model_fit_mmcp = model_fit_mmcp[which.min(model_fit_mmcp[,3]), , drop = FALSE]
  t1 = Sys.time()
  
  # NSD (nonlinear least squares) model
  model_fit_nsd = fit_NSD_migration(data_for_migration[, c(1,2,7)], units = "km")
  t2 = Sys.time()
  
  # FPT (penalized contrast) model
  model_fit_fpt = fit_fpt_penalized_contrast(data_for_migration[, c(1,2,7)], seq(5, 100, 5))
  t3 = Sys.time()
  
  # Piecewise regression (Bayesian) model
  model_fit_pwr = fit_NSD_piecewise(data_for_migration)
  t4 = Sys.time()
  
  # Mechanistic range shift analysis
  model_fit_mrsa = fit_marcher(data_for_migration)
  t5 = Sys.time()
  
  # Behavioral change point analysis
  model_fit_bcpa = fit_bcpa(data_for_migration[, c(1,2,7)]) %>% sort
  t6 = Sys.time()
  
  # Bayesian partitioning of Markov models
  model_fit_bpmm = fit_bpmm(data_for_migration[, c(1, 2, 7)])
  t7 = Sys.time()
  
  comp_times = as.numeric(c(t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6))
  this_result = c(model_fit_mmcp, model_fit_nsd, model_fit_fpt, model_fit_pwr, model_fit_mrsa, model_fit_bcpa, model_fit_bpmm, comp_times)
  all_results = rbind(all_results, this_result)
  
  message("Completed simulation ", i)
}