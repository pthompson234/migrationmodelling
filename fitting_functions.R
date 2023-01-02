# Simulate a migratory movement
#
# n_timesteps: length of the simulated track (NOTE: if percent_removed > 0, this is NOT necessarily the number of total data points!)
# migration_starts: vector of start times for the switch to migratory behavior
# migration_ends: vector of end times for migratory behavior (will switch back to sedentary behavior)
# percent_removed: number between 0 and 1; how many of the timesteps should we randomly remove (simulates missing data)
# sl_migration: mean step length on migration
# sl_resident: mean step length during "sedentary" movement
# kappa_migration: directional autocorrelation parameter on migration (turn angle distribution is uniform for sedentary movement)
# kappa_resident: directional autocorrelation parameter outside of migration
# return_code2_only: do we only return code 2's (see function body for description) or return all non-missing data points
# Returns a data.frame ready to be passed to the fit_migration_TMB function (see documentation for that function)
sim_migratory_movement = function(n_timesteps, 
                                  migration_starts, 
                                  migration_ends, 
                                  percent_removed = 0.5,
                                  sl_migration = 5, 
                                  sl_resident = 1,
                                  kappa_migration = 10,
                                  kappa_resident = 0,
                                  return_code2_only = TRUE) {
  
  # Each step gets a "Code" representing its applicability to model fitting (i.e., does it include any missing data)
  # -2: "missing" point that is included in data frame temporarily
  # -1: "available" point simulated for use-available modelling
  # 0: endpoint of a step with no begin point
  # 1: endpoint of a step with only one previous consecutive point
  # 2: endpoint of a "full" step with two consecutive previous points
  # 3: first point of the dataset
  # 4: endpoint of step beginning at the fisrt point of the dataset
  
  require(circular) # for von Mises distribution
  
  all_df = data.frame(t = 0, x = 0, y = 0, sl = 0, he = 0, ta = 0, Code = 3)
  
  for (ts in 1:n_timesteps) {
    is_it_migration = (sum(ts > migration_starts) + sum(ts > migration_ends)) %% 2
    if (is_it_migration) {
      # Simulate a migratory step
      r_steplength = rexp(1, 1 / sl_migration)
      r_turnangle = as.numeric(rvonmises(1, 0, kappa_migration))
    } else {
      # Simulate a "sedentary" step
      r_steplength = rexp(1, 1 / sl_resident)
      r_turnangle = as.numeric(rvonmises(1, 0, kappa_resident))
    }
    
    r_heading = (all_df$he[ts] + r_turnangle) %% (2*pi)
    r_x = r_steplength * cos(r_heading) + all_df$x[ts]
    r_y = r_steplength * sin(r_heading) + all_df$y[ts]
    
    if (runif(1) < percent_removed) {
      # include this point with Code == -2
      r_code = -2
    } else if (all_df$Code[ts] == 0) {
      r_code = 1
    } else if (all_df$Code[ts] == 3) {
      # this is the second step
      r_code = 4
    } else if (all_df$Code[ts] > 0) {
      # if code is equal to 1, 2, or 4, then we have a good step
      r_code = 2
    } else {
      r_code = 0
    }
    
    all_df = rbind(all_df, cbind(t = ts, x = r_x, y = r_y, sl = r_steplength, he = r_heading, ta = r_turnangle, Code = r_code))
    
  }
  
  revised_df = all_df[3:nrow(all_df), ] # get only points where we can get (x2, y2)
  revised_df$x1 = all_df$x[2:(nrow(all_df) - 1)]
  revised_df$y1 = all_df$y[2:(nrow(all_df) - 1)]
  revised_df$x2 = all_df$x[1:(nrow(all_df) - 2)]
  revised_df$y2 = all_df$y[1:(nrow(all_df) - 2)]
  
  if (!return_code2_only) return(revised_df[revised_df$Code >= 0, c("x", "y", "x1", "y1", "x2", "y2", "t", "Code")])
  revised_df[revised_df$Code == 2, c("x", "y", "x1", "y1", "x2", "y2", "t")]
  
}

# Fit the one-migration model (with step function and fixed time) using TMB
#
# data_in: matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# times: 2-vector (should be ordered) of start times (t1, t2)
# return_LL: do we return LL_2r to the global environment
fit_migration_c1_TMB = function(data_in, times, return_LL = FALSE) {
  
  # require(optimx) # don't think I need for now but could come back later
  require(TMB)
  
  Y_steps = data_in[,c(1,2)] - data_in[,c(3,4)]
  steplengths = sqrt(rowSums(Y_steps^2))
  avg_speed = mean(steplengths)
  max_speed = max(steplengths)
  Y_prev = data_in[,c(3,4)] - data_in[,c(5,6)]
  headings = atan2(Y_steps[, 2], Y_steps[, 1])
  prev_headings = atan2(Y_prev[,2], Y_prev[,1])
  turn_angles = headings - prev_headings
  
  filename = "migration_c1"
  compile(paste0(filename, ".cpp"))
  
  data_model = list(steplengths = steplengths, turn_angles = turn_angles, time_indices = data_in[, 7], t1 = times[1], t2 = times[2])
  pars_model = list(r0 = 1, r1 = 1, k0 = 1, k1 = 1) # these values actually don't matter
  dyn.load(dynlib(filename))
  
  LL_1rf = MakeADFun(data = data_model, parameters = pars_model, DLL = filename)
  
  model_fit = nlminb(LL_1rf$par, LL_1rf$fn, LL_1rf$gr, lower = rep(0, 4), upper = c(max_speed, max_speed, 1e5, 1e5), control = list(eval.max = 1000, iter.max = 1000))
  
  ses = sdreport(LL_1rf)$sd
  
  if (return_LL) LL_1rf <<- LL_1rf
  
  data.frame(model = "1rf", t1 = times[1], t2 = times[2], parname = names(pars_model), NLL = model_fit$objective, AIC = 2*model_fit$objective + 2*length(model_fit$par), BIC = 2*model_fit$objective + log(nrow(data_in))*length(model_fit$par), estimate = model_fit$par, std_err = ses, CL = model_fit$par - 1.96*ses, CU = model_fit$par + 1.96*ses, convergence = model_fit$message)
  
}

# Fit the two-migration model (with step function and fixed time and only one set of migratory parameters) using TMB
#
# data_in: matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# times: 4-vector (should be ordered) of start times (t1, t2, t3, t4)
# return_LL: do we return LL_2r to the global environment
fit_migration_c2_TMB = function(data_in, times, return_LL = FALSE) {
  
  # require(optimx) # don't think I need for now but could come back later
  require(TMB)
  
  Y_steps = data_in[,c(1,2)] - data_in[,c(3,4)]
  steplengths = sqrt(rowSums(Y_steps^2))
  avg_speed = mean(steplengths)
  max_speed = max(steplengths)
  Y_prev = data_in[,c(3,4)] - data_in[,c(5,6)]
  headings = atan2(Y_steps[, 2], Y_steps[, 1])
  prev_headings = atan2(Y_prev[,2], Y_prev[,1])
  turn_angles = headings - prev_headings
  
  filename = "migration_c2"
  compile(paste0(filename, ".cpp"))
  
  data_model = list(steplengths = steplengths, turn_angles = turn_angles, time_indices = data_in[, 7], t1 = times[1], t2 = times[2], t3 = times[3], t4 = times[4])
  pars_model = list(r0 = 1, r1 = 1, k0 = 1, k1 = 1) # these values actually don't matter
  dyn.load(dynlib(filename))
  
  LL_2rf1p = MakeADFun(data = data_model, parameters = pars_model, DLL = filename)
  
  model_fit = nlminb(LL_2rf1p$par, LL_2rf1p$fn, LL_2rf1p$gr, lower = rep(0,4), upper = c(max_speed, max_speed, 1e5, 1e5), control = list(eval.max = 1000, iter.max = 1000))
  
  ses = sdreport(LL_2rf1p)$sd
  
  if (return_LL) LL_2rf1p <<- LL_2rf1p
  
  data.frame(model = "2rf1p", t1 = times[1], t2 = times[2], t3 = times[3], t4 = times[4], parname = names(pars_model), NLL = model_fit$objective, AIC = 2*model_fit$objective + 2*length(model_fit$par), BIC = 2*model_fit$objective + log(nrow(data_in))*length(model_fit$par), estimate = model_fit$par, std_err = ses, CL = model_fit$par - 1.96*ses, CU = model_fit$par + 1.96*ses, convergence = model_fit$message)
  
}

# Make a matrix of possible migration times based on a set of inputs
#
# starts = vector of start times
# ends = vector of end times
# starts2 = vector of start times for second migration
# ends2 = vector of end times for second migration
# n_migrations: how many migrations are we getting
# min_diffs: minimum required difference between any start-end (or start2-end2) pair. Either a number (applies to all) or a vector (t2-t1, t3-t2, t4-t3)
#
# returns a matrix with 2 * n_migrations columns [start1, end1, start2, end2] or [start1, end1]
make_times_matrix = function(starts, ends = starts, starts2 = ends, ends2 = starts2, n_migrations = 1, min_diff = 0) {
  
  if (n_migrations == 1) {
    matrix_all = cbind(rep(starts, times = length(ends)), rep(ends, each = length(starts)))
    
    return(matrix_all[matrix_all[,2]-matrix_all[,1] > min_diff, ])
  }
  
  if (length(min_diff) == 1) min_diff = rep(min_diff, 3) # turn min_diff into a vector of appropriate length if it isn't already one
  
  single_times1 = make_times_matrix(starts, ends, n_migrations = 1)
  single_times2 = make_times_matrix(starts2, ends2, n_migrations = 1)
  matrix_all = cbind(apply(single_times1, 2, rep, times = nrow(single_times2)),
                     apply(single_times2, 2, rep, each = nrow(single_times1)))
  matrix_all = matrix_all[((matrix_all[,2] - matrix_all[,1]) > min_diff[1]) & ((matrix_all[,3] - matrix_all[,2]) > min_diff[2]) & ((matrix_all[,4] - matrix_all[,3]) > min_diff[3]), ]
  return(matrix_all)
  
}

# Fits the multivariate Gaussian migration model using TMB (r(t) is a step function with fixed bounds)
#
# times_list: matrix or data.frame with 2 or 4 columns (one for start time, one for end time)
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# ...: additional arguments to fitting function
#
# Returns a model fit 
estimate_times = function(times_list, data_in, ...) {
  
  # Remove items from times_list that are greater than max(data_in[,7]) (i.e., one of the times is too large)
  times_list = times_list[times_list[, ncol(times_list)] <= max(data_in[, 7]), ]
  
  result = data.frame()
  
  for (i in 1:nrow(times_list)) {
    if (ncol(times_list) == 2) {
      fit = tryCatch(fit_migration_c1_TMB(data_in, as.numeric(times_list[i, ]), ...),
                     error = function(e) {
                       data.frame(NLL = Inf, estimate = rep(NA, 4), convergence = as.character(e))
                     })
      
      result = rbind(result, cbind(t1 = times_list[i, 1], t2 = times_list[i, 2], NLL = fit$NLL[1], r0 = fit$estimate[1], 
                                   r1 = fit$estimate[2], k0 = fit$estimate[3], k1 = fit$estimate[4], msg = fit$convergence[1]))
    } else {
      fit = tryCatch(fit_migration_c2_TMB(data_in, as.numeric(times_list[i, ]), ...),
                     error = function(e) {
                       data.frame(NLL = Inf, estimate = rep(NA, 4), convergence = as.character(e))
                     })
      
      result = rbind(result, cbind(t1 = times_list[i, 1], t2 = times_list[i, 2], t3 = times_list[i, 3], t4 = times_list[i, 4], 
                                   NLL = fit$NLL[1], r0 = fit$estimate[1], r1 = fit$estimate[2], k0 = fit$estimate[3], 
                                   k1 = fit$estimate[4], msg = fit$convergence[1]))
    }
  }
  
  result
  
}

# Run the estimate_times function twice with a coarse (times1) and fine (times2) temporal resolution to more efficiently identify the optimum
#
# times1: matrix or data.frame with 2 or 4 columns (one for start time, one for end time)
# interval2: integer vector of intervals to test (the length of the vector is the # of grid levels minus 1)
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# NLL_tol: how close does the NLL from rows of times1 have to be to the optimum to be included in the subset of times2?
# NLL_include: how many of the best NLL's to include for the next step? If not NULL, overrides NLL_tol.
# buffers: numbers - see line with big make_times_matrix call
# bounds: minimum and maximum temporal extents (2-vector; any additional components will be ignored)
# min_diff: for make_times_matrix (see documentation there)
# ...: additional arguments to fitting function
#
# Returns a data frame similar to estimate_times
estimate_times_nlevels = function(times1, intervals = NULL, data_in = NULL, NLL_tol = 0, NLL_include = NULL, buffers = c(1, 1, 1, 1), bounds = c(-Inf, Inf), min_diff = 0, ...) {
  
  if (is.null(data_in)) stop("need data")
  if (is.null(intervals)) return(estimate_times(times1, data_in, ...))
  
  if (!is.null(NLL_include) && NLL_include < 0.5) stop("Invalid 'NLL_include' value. Should be an integer greater than 1.")
  
  estimates = estimate_times(times1, data_in, ...)
  if (NLL_include > nrow(estimates)) NLL_include = nrow(estimates) # so we don't have more rows than we need!
    
  if (any(order(intervals) != (length(intervals):1))) stop("intervals not ordered properly")
  
  seq_fun = function(x, ind_1 = 2, ind_2 = 1) {
    seq(x - intervals[ind_1 - 1] * buffers[ind_2], x + intervals[ind_1 - 1] * buffers[ind_2], intervals[ind_1])
  }
  
  for (int in 2:length(intervals)) {
    if (is.null(NLL_include)) {
      result_refined = estimates[estimates$NLL <= (min(as.numeric(estimates$NLL)) + NLL_tol), ]
    } else {
      result_refined = estimates[order(as.numeric(estimates$NLL))[1:round(NLL_include)], ]
    }
    
    t1_seq = unique(as.numeric(sapply(X = as.numeric(result_refined$t1), FUN = seq_fun, ind_1 = int, ind_2 = 1)))
    t2_seq = unique(as.numeric(sapply(X = as.numeric(result_refined$t2), FUN = seq_fun, ind_1 = int, ind_2 = 2)))
    
    if (ncol(times1) == 2) {
      # only estimating one migration
      times2 = make_times_matrix(starts = t1_seq[order(t1_seq)], ends = t2_seq[order(t2_seq)], n_migrations = 1, min_diff = min_diff)
    } else {
      # estimating two migrations
      t3_seq = unique(as.numeric(sapply(X = as.numeric(result_refined$t3), FUN = seq_fun, ind_1 = int, ind_2 = 3)))
      t4_seq = unique(as.numeric(sapply(X = as.numeric(result_refined$t4), FUN = seq_fun, ind_1 = int, ind_2 = 4)))
      
      times2 = make_times_matrix(starts = t1_seq[order(t1_seq)], ends = t2_seq[order(t2_seq)], starts2 = t3_seq[order(t3_seq)], 
                                 ends2 = t4_seq[order(t4_seq)], n_migrations = 2, min_diff = min_diff)
    }
    
    # remove times in times2 that are outside of bounds
    times2 = times2[(apply(times2, 1, min) >= bounds[1]) & (apply(times2, 1, max) <= bounds[2]), ]
    
    estimates = estimate_times(times2, data_in, ...)
  }
  
  estimates
  
}

# Calculate net squared displacement (NSD) from data
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y]. Note that full steps and time indices are not needed here!
#
# Returns a vector of NSD values for each row in data_in
net_squared_displacement = function(data_in) {
  sapply(X = 1:nrow(data_in), FUN = function(row) {
    (data_in[row, 1] - data_in[1, 1])^2 + (data_in[row, 2] - data_in[1, 2])^2
  })
}

# Estimate the beginning and end of migration using the Net Squared Displacement (NSD) approach from Bunnefeld et al. (2011)
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, difftime/# of days from start]. Note that full steps are not needed here!
# units: what units are the data in? Right now code handles "m" (meters) and "km" (kilmeters). The only reason this matters is because if units == "m" the NSD values become ridiculously large and it's problematic. But, if units == "km" we don't want to divide unnecessarily.
# Output: 
fit_NSD_migration = function(data_in, 
                             units = "m") {
  # Note we use 1000^2 here because NSD is a squared unit
  factor_div = ifelse(units == "m", 1000^2, 1)
  NSD = net_squared_displacement(data_in) / factor_div
  t = data_in[, 3]
  
  model_nsd = nls(NSD ~ delta / (1 + exp((theta - t)/phi)), start = list(delta = median(NSD), theta = median(t), phi = 10), trace = 0, lower = list(delta = 0, theta = min(t), phi = 1), upper = list(delta = max(NSD), theta = max(t), phi = max(t)-min(t)), algorithm = "port")
  
  # Get model parameters and convert them into t1, t2, and intensity. For now, leave out CI's because I don't think you can extrapolate CI's for t1 and t2 in the way that I'm imagining, because the parameters are not independent.
  pars_orig = coef(model_nsd)
  # ses_orig = summary(model)$coefficients[,2]
  # CI_lower = pars_orig - 1.96 * ses_orig
  # CI_upper = pars_orig + 1.96 * ses_orig
  # pars_revised = c(t1 = pars_orig[2] - 2 * pars_orig[3], t2 = pars_orig[2] + 2 * pars_orig[3], intensity = pars_orig[1])
  # CI_lower_revised = c(start = CI_lower[2] - 2 * CI_upper[3], end = CI_lower[2] + 2 * CI_lower[3], intensity = CI_lower[1])
  # CI_upper_revised = c(start = CI_upper[2] - 2 * CI_lower[3], end = CI_upper[2] + 2 * CI_upper[3], intensity = CI_upper[1])
  
  c(pars_orig, t1 = pars_orig[2] - 2 * pars_orig[3], t2 = pars_orig[2] + 2 * pars_orig[3])
}

# Helper function to get the negative log-likelihood for the step-length based calving model (Model 1 from DeMars et al., 2013) only for the period after calving (t > t1)
#
# step_lengths: array of step lengths
# times: array of times corresponding to the step lengths
# t1: beginning of calving period
# k: duration of calving period
nll_approx_k = function(step_lengths, times, t1, k) {
  
  # manually calculate mean step length for non-calving period (it will optimize)
  lambda = mean(step_lengths[times <= t1 | times > t1 + k])
  
  # get individual mean step length values
  rho_values = (times - t1) * lambda / k
  rho_values[rho_values > lambda | rho_values <= 0] = lambda
  # rho_values[rho_values <= lambda] = 1 # Only uncomment this if we want to test another kind of model
  data.frame(NLL = -sum(dexp(step_lengths, 1 / rho_values, log = TRUE)), t1 = t1, k = k, lambda = lambda)
  
}

# Get optimal timing for calving period
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# times_array: 2-column matrix or data.frame with values for t1 and t2 (t1 + k)
fit_calving = function(data_in, times_array) {
  
  step_lengths = sqrt((data_in[,1] - data_in[,3])^2 + (data_in[,2] - data_in[,4])^2)
  times = data_in[,7]
  
  # as per DeMars et al. (2013) cut off top percentile of step lengths to fix the mean problem
  top_perc_sl = quantile(step_lengths, 0.99)
  include_step = step_lengths < top_perc_sl
  step_lengths = step_lengths[include_step]
  times = times[include_step]
  
  # Remove instances of times from times_array if those times aren't recorded in the data
  times_array = times_array[(times_array[,1] %in% times) & (times_array[,2] %in% times), ]
  
  t1_array = times_array[,1]
  k_array = times_array[,2] - times_array[,1]
  
  model_fits = sapply(1:nrow(times_array), FUN = function(row) {
    nll_approx_k(step_lengths, times, t1_array[row], k_array[row])
  })
  as.data.frame(t(model_fits))
  
}

# Fit the model from Le Corre et al. (2014) to migration data
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, difftime/# of days from start]. Note that full steps are not needed here!
# fpt_radii: what radii do we want to test for the first-passage-time analysis? Note that the final analysis will only use one of these radii
#
# Returns a 2-vector of t1 and t2
fit_fpt_penalized_constant = function(data_in,
                                      fpt_radii) {
  require(adehabitatLT)
  
  # Make data into ltraj format for FPT analysis. Note we have to turn the times into dates but it doesn't really matter what they are
  data_ltraj = as.ltraj(xy = data_in[,1:2], date = as.POSIXct("2000-01-01") + data_in[,3] * 60 * 60 * 24, id = "data")
  data_fpt_profile = fpt(data_ltraj, fpt_radii, units = "days")
  # Get the variance of log(FPT) for each radius and the maximum value of this statistic will become the radius we use
  data_var_fpt = varlogfpt(data_fpt_profile, graph = FALSE)
  chosen_fpt_radius_index = which.max(as.numeric(data_var_fpt[1, ]))
  
  # Get FPT time series for the actual data
  fpt_time_series = data_fpt_profile[[1]][,chosen_fpt_radius_index]
  # This time series will have many NA values so we change them to the maximum time
  fpt_time_series[is.na(fpt_time_series)] = max(data_in[,3])
  
  # Use the penalized constant method from Lavielle (2005) to identify the optimal breakpoints
  penalized_constants = lavielle(fpt_time_series, Lmin = 1, Kmax = 3, type = "mean")
  optimal_breakpoints = findpath(penalized_constants, K = 3, plotit = FALSE) %>% unlist
  optimal_times = data_in[optimal_breakpoints, 3]
  
  # Return times corresponding to the breakpoints in the path analysis
  c(t1 = mean(optimal_times[2:3]), t2 = mean(optimal_times[4:5]))
}

# Fit the model from Madon & Hingrat (2014) to migration data
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start, Code]. Can include non-code 2's although anything 0 or lower will be removed.
# fpt_radii: what radii do we want to test for the first-passage-time analysis? Note that the final analysis will only use one of these radii
#
# Returns a 2-vector of t1 and t2
fit_dailydist_pelt = function(data_in) {
  # Reduce the data to steps with beginning points only
  data_in = data_in[data_in$Code > 0 & data_in$Code != 3, ]
  
  distances = sqrt((data_in[,2] - data_in[,4])^2 + (data_in[,1] - data_in[,3])^2)
  # Copying the methods from Madon & Hingrat to my best ability
  changepoints_pelt = cpt.var(distances, penalty = "SIC", method = "PELT", Q = 2)
  changepoint_values = changepoints_pelt@cpts[1:2] # only take the first two; the third point is the ending value if it's there
  c(t1 = data_in[changepoint_values[1], 7], t2 = data_in[changepoint_values[2], 7])
}

# Given a dataset with a known optimum, calculate bootstrapped confidence intervals
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# optim_pars: named vector (important: the vector must be named - this helps determine what model we are fitting)
#
# Output: a data.frame just like data_in (or at least, with the necessary columns in indices 1-7)
sim_sl_ta = function(data_in, pars) {
  
  require(circular) # von Mises distribution
  
  rho_pars = rep(pars['r0'], nrow(data_in))
  kappa_pars = rep(pars['k0'], nrow(data_in))
  
  rho_pars[data_in[,7] > pars['t1'] & data_in[,7] < pars['t2']] = pars['r0'] + pars['r1']
  kappa_pars[data_in[,7] > pars['t1'] & data_in[,7] < pars['t2']] = pars['k0'] + pars['k1']
  
  # update further if we have 2 migrations
  if (!is.na(pars['t3']) & is.na(pars['r2'])) {
    rho_pars[data_in[,7] > pars['t3'] & data_in[,7] < pars['t4']] = pars['r0'] + pars['r1']
    kappa_pars[data_in[,7] > pars['t3'] & data_in[,7] < pars['t4']] = pars['k0'] + pars['k1']
  } else if (!is.na(pars['r2'])) {
    rho_pars[data_in[,7] > pars['t3'] & data_in[,7] < pars['t4']] = pars['r0'] + pars['r2']
    kappa_pars[data_in[,7] > pars['t3'] & data_in[,7] < pars['t4']] = pars['k0'] + pars['k2']
  }
  
  # get step lengths and turning angles
  r_steplengths = rexp(nrow(data_in), 1 / rho_pars)
  r_turn_angles = sapply(kappa_pars, function(k) {
    as.numeric(rvonmises(1, 0, k))
  })
  
  # get actual headings, which we add turning angles on to to get the new headings
  prev_headings = atan2(data_in[,4] - data_in[,6], data_in[,3] - data_in[,5])
  
  r_headings = prev_headings + r_turn_angles
  r_x = data_in[,3] + r_steplengths * cos(r_headings)
  r_y = data_in[,4] + r_steplengths * sin(r_headings)
  
  simmed_data = data_in
  simmed_data[,1] = r_x
  simmed_data[,2] = r_y
  
  simmed_data
  
}

# Given a dataset with a known optimum, calculate bootstrapped confidence intervals
#
# data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
# optim_pars: named vector (important: the vector must be named - this helps determine what model we are fitting)
# times_matrix: for estimate_times_nlevels, times matrix used for model fitting (should be the same one as used to obtain optim_pars from the original fit)
# n_iter: how many simulations to run (how many bootstraps?)
# parallel: how many cores to distribute the job over? If 1 or less, we don't actually use the parallel implementation
# ...: additional arguments to esitmate_times_nlevels
#
# Output: a matrix of confidence bounds for each of the values in optim_pars
bootstrap_CI = function(data_in, optim_pars, times_matrix, n_iter = 100, parallel = 1, confidence_level = 0.95, ...) {
  
  require(circular)
  
  if (parallel > 1) {
    require(doParallel)
    registerDoParallel(cores = parallel)
    
    all_sim_results = foreach(i = 1:n_iter, .combine = rbind, .errorhandling = "stop") %dopar% {
      simmed_steps = sim_sl_ta(data_in, optim_pars)
      
      model_fit = estimate_times_nlevels(times1 = times_matrix, data_in = simmed_steps, ...)
      
      model_fit[which.min(model_fit$NLL), ]
    }
  } else {
    all_sim_results = data.frame()
    
    for (i in 1:n_iter) {
      simmed_steps = sim_sl_ta(data_in, optim_pars)
      model_fit = estimate_times_nlevels(times1 = times_matrix, data_in = simmed_steps, ...)
      all_sim_results = rbind(all_sim_results, model_fit[which.min(model_fit$NLL), ])
    }
  }
  
  CIs_matrix = matrix(ncol = 2, nrow = length(optim_pars), 0)
  for (row in 1:length(optim_pars)) {
    estimates = as.numeric(all_sim_results[, names(optim_pars)[row]])
    CIs_matrix[row, 1] = quantile(estimates, (1 - confidence_level) / 2)
    CIs_matrix[row, 2] = quantile(estimates, 1 - (1 - confidence_level) / 2)
  }
  
  CIs_matrix
}