library(adehabitatLT)
library(changepoint)
library(bcpa)
library(tidyverse)
library(mcp)
library(marcher)
library(rjags)
library(Rcpp)

sourceCpp("C:/Users/pthom/Documents/School/THESIS/Ch4_5_migration/code/Github/cpp_functions.cpp")
source("C:/Users/pthom/Documents/School/THESIS/Ch4_5_migration/code/Github/fitting_functions.R")

N_SIM = 55 # How many random migrations to simulate? In reality we want 50; this is coverage for the occasional failed convergence from the NSD model.

# True values for all model parameters
t1_MIG = 70
t2_MIG = 100

N_DAYS = 200
N_STEPS_DESIRED = 150
GPS_ERROR = 0.5

RHO_0 = 5
RHO_1 = 40
KAPPA_0 = 0
KAPPA_1 = 0.5

BUFFER = 7
FPT_RADII = seq(5, 100, 5)

for (i in 1:N_SIM) {
  # Simulate random data
  data_for_migration = sim_migratory_movement(N_DAYS, t1_MIG, t2_MIG, percent_removed = (1 - N_STEPS_DESIRED / N_DAYS) / 2, sl_migration = RHO_0 + RHO_1, sl_resident = RHO_0, kappa_migration = KAPPA_1, kappa_resident = KAPPA_0, error = GPS_ERROR, return_code2_only = FALSE)
  
  data_code2only = data_for_migration[data_for_migration$Code == 2, -8] %>% as.matrix # for our model
  
  # Fit our model
  times_matrix = make_times_matrix(c(seq(min(data_code2only[, 7]), max(data_code2only[, 7]), by = 14), max(data_code2only[,7])), n_migrations = 1, min_diff = BUFFER)
  
  message("Starting model fit for simulation ", i)
  
  t1 = Sys.time()
  
  model_fit_thompson = estimate_times_nlevels(times_matrix, c(14, 7, 3, 1), data_code2only, mod_type = "r1k1", NLL_include = 5, cpp = TRUE, bounds = c(min(data_code2only[, 7]), max(data_code2only[, 7])), min_diff = BUFFER)
  model_fit_thompson = model_fit_thompson[which.min(model_fit_thompson[,3]), , drop = FALSE]
  print("Done")
  
  t2 = Sys.time()
  
  # Fit the NSD model from Bunnefeld et al (2011)
  model_fit_bunnefeld = tryCatch(fit_NSD_migration(data_for_migration[, c(1,2,7)], units = "km"), error = function(e) {
    return(c(delta = NA, theta = NA, phi = NA, t1.theta = NA, t2.theta = NA))
  })
  
  t3 = Sys.time()
  
  # Fit the FPT model from Le Corre et al. (2014)
  model_fit_lecorre = fit_fpt_penalized_contrast(data_for_migration[, c(1,2,7)], FPT_RADII)
  
  t4 = Sys.time()
  
  # Fit the PELT model using daily distances from Madon & Hingrat (2014)
  model_fit_madon = fit_dailydist_pelt(data_for_migration)
  
  t5 = Sys.time()
  
  # Fit the Bayesian piecewise regression model from Wolfson et al. (2022)
  model_fit_wolfson = fit_NSD_piecewise(data_for_migration)
  
  t6 = Sys.time()
  
  # Fit the marcher model (Gurarie et al., 2017)
  model_fit_marcher = fit_marcher(data_for_migration)
  
  t7 = Sys.time()
  
  diff_times = data.frame(dt1 = as.numeric(difftime(t2, t1, units = "secs")),
                          dt2 = as.numeric(difftime(t3, t2, units = "secs")),
                          dt3 = as.numeric(difftime(t4, t3, units = "secs")),
                          dt4 = as.numeric(difftime(t5, t4, units = "secs")),
                          dt5 = as.numeric(difftime(t6, t5, units = "secs")),
                          dt6 = as.numeric(difftime(t7, t6, units = "secs")))
  
  result = cbind(N = i, model_fit_thompson, t(model_fit_bunnefeld), t(model_fit_lecorre), t(model_fit_madon), t(model_fit_wolfson), t(model_fit_marcher), diff_times)
  all_results = rbind(all_results, result)
  
  message("Completed simulation ", i)
}