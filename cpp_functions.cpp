#include <Rcpp.h>
using namespace Rcpp;

// for solving MLE for kappa; approximates the ratio of inverse Bessel functions needed to acquire MLE for kappa. Only for internal use.
//
// x: value representing mean cosine of angles; should be between 0 and 1 as a result (negative values indicate a kappa MLE of 0)
//
// Returns a non-negative real number representing the MLE for kappa
double A1Inv(double x) {
  if (x >= 0 & x < 0.53) {
    return 2 * x + x*x*x + (5 * x*x*x*x*x)/6;
  } else if (x < 0.85) {
    return -0.4 + 1.39 * x + 0.43/(1 - x);
  } else {
    return 1/(x*x*x - 4 * x*x + 3 * x);
  }
}

// obtain the (approximate) analytical MLE for the kappa parameter of the von Mises distribution
//
// angles: vector of angular values (can be any real number)
//
// Returns a non-negative real integer representing the analytical MLE for kappa
double mle_vm_kappa(NumericVector angles) {
  int n = angles.size();
  
  double v = mean(cos(angles));
  double kappa_est = 0;
  if (v > 0) {
    kappa_est = A1Inv(v);
  }
  return kappa_est;
}

// get the angle generated from x and y by calculating atan(y/x) with some corrections
//
// x vector of x coordinates
// y vector of y coordinates
//
// returns a vector of angles (in radians; [-pi/2, 3pi/2])
NumericVector atan2_cpp(NumericVector x, NumericVector y) {
  NumericVector result = ifelse(x < 0, atan(y / x) + M_PI, atan(y / x));
  result[x == 0 & y > 0] = M_PI_2;
  result[x == 0 & y < 0] = -1 * M_PI_2;
  result[x == 0 & y == 0] = 0.0;
  return result;
}

// obtain MLE's (and negative log likelihood) for the c = 1 migration model
//
// data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
// times: vector of length 2c representing the (fixed) start and end time of migration
//
// Returns a vector of the following form: NULL, rho0, rho1, kappa0, kappa1
// [[Rcpp::export]]
NumericVector fit_1migration(NumericMatrix data_in, NumericVector times) {
  
  NumericVector x_steps = data_in( _ , 0) - data_in( _ , 2);
  NumericVector y_steps = data_in( _ , 1) - data_in( _ , 3);
  NumericVector steplengths = sqrt(x_steps * x_steps + y_steps * y_steps);
  NumericVector headings = atan2_cpp(x_steps, y_steps);
  NumericVector x_prev_steps = data_in( _ , 2) - data_in( _ , 4);
  NumericVector y_prev_steps = data_in( _ , 3) - data_in( _ , 5);
  NumericVector prev_headings = atan2_cpp(x_prev_steps, y_prev_steps);
  NumericVector turn_angles = headings - prev_headings;
  NumericVector t_values = data_in( _ , 6);
  
  NumericVector steplengths_0 = steplengths[t_values <= times(0) | t_values > times(1)];
  NumericVector steplengths_1 = steplengths[t_values > times(0) & t_values <= times(1)];
  NumericVector turnangles_0 = turn_angles[t_values <= times(0) | t_values > times(1)];
  NumericVector turnangles_1 = turn_angles[t_values > times(0) & t_values <= times(1)];
  
  double rho_0 = mean(steplengths_0);
  double rho_1 = mean(steplengths_1) - rho_0;
  if (rho_1 < 0) {
    rho_1 = 0;
  }
  double kappa_0 = mle_vm_kappa(turnangles_0);
  double kappa_1 = mle_vm_kappa(turnangles_1) - kappa_0;
  if (kappa_1 < 0) {
    kappa_1 = 0;
  }
  
  double NLL_rho0 = sum(steplengths_0 / rho_0 - log(1 / rho_0));
  double NLL_rho1 = sum(steplengths_1 / (rho_0 + rho_1) - log(1 / (rho_0 + rho_1)));
  double NLL_kappa0 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0, 0, 1)) - kappa_0 * cos(turnangles_0));
  double NLL_kappa1 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0 + kappa_1, 0, 1)) - (kappa_0 + kappa_1) * cos(turnangles_1));
  
  return {NLL_rho0 + NLL_rho1 + NLL_kappa0 + NLL_kappa1, rho_0, rho_1, kappa_0, kappa_1};
  
}

// obtain MLE's (and negative log likelihood) for the c = 2 migration model
//
// data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
// times: vector of length 2c representing the (fixed) start and end time of migration
//
// Returns a vector of the following form: NULL, rho0, rho1, kappa0, kappa1
// [[Rcpp::export]]
NumericVector fit_2migrationr1k1(NumericMatrix data_in, NumericVector times) {

  NumericVector x_steps = data_in( _ , 0) - data_in( _ , 2);
  NumericVector y_steps = data_in( _ , 1) - data_in( _ , 3);
  NumericVector steplengths = sqrt(x_steps * x_steps + y_steps * y_steps);
  NumericVector headings = atan2_cpp(x_steps, y_steps);
  NumericVector x_prev_steps = data_in( _ , 2) - data_in( _ , 4);
  NumericVector y_prev_steps = data_in( _ , 3) - data_in( _ , 5);
  NumericVector prev_headings = atan2_cpp(x_prev_steps, y_prev_steps);
  NumericVector turn_angles = headings - prev_headings;
  NumericVector t_values = data_in( _ , 6);
  
  NumericVector steplengths_0 = steplengths[t_values <= times(0) | t_values > times(3) | (t_values > times(1) & t_values <= times(2))];
  NumericVector steplengths_1 = steplengths[(t_values > times(0) & t_values <= times(1)) | (t_values > times(2) & t_values <= times(3))];
  NumericVector turnangles_0 = turn_angles[t_values <= times(0) | t_values > times(3) | (t_values > times(1) & t_values <= times(2))];
  NumericVector turnangles_1 = turn_angles[(t_values > times(0) & t_values <= times(1)) | (t_values > times(2) & t_values <= times(3))];
  
  double rho_0 = mean(steplengths_0);
  double rho_1 = mean(steplengths_1) - rho_0;
  if (rho_1 < 0) {
    rho_1 = 0;
  }
  double kappa_0 = mle_vm_kappa(turnangles_0);
  double kappa_1 = mle_vm_kappa(turnangles_1) - kappa_0;
  if (kappa_1 < 0) {
    kappa_1 = 0;
  }
  
  double NLL_rho0 = sum(steplengths_0 / rho_0 - log(1 / rho_0));
  double NLL_rho1 = sum(steplengths_1 / (rho_0 + rho_1) - log(1 / (rho_0 + rho_1)));
  double NLL_kappa0 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0, 0, 1)) - kappa_0 * cos(turnangles_0));
  double NLL_kappa1 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0 + kappa_1, 0, 1)) - (kappa_0 + kappa_1) * cos(turnangles_1));
  
  return {NLL_rho0 + NLL_rho1 + NLL_kappa0 + NLL_kappa1, rho_0, rho_1, kappa_0, kappa_1};
  
}

// obtain MLE's (and negative log likelihood) for the c = 2 migration model with unique kappa's and rho's for each migration
//
// data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
// times: vector of length 2c representing the (fixed) start and end time of migration
//
// Returns a vector of the following form: NULL, rho0, rho1, rho2, kappa0, kappa1, kappa2
// [[Rcpp::export]]
NumericVector fit_2migrationrCkC(NumericMatrix data_in, NumericVector times) {
  
  NumericVector x_steps = data_in( _ , 0) - data_in( _ , 2);
  NumericVector y_steps = data_in( _ , 1) - data_in( _ , 3);
  NumericVector steplengths = sqrt(x_steps * x_steps + y_steps * y_steps);
  NumericVector headings = atan2_cpp(x_steps, y_steps);
  NumericVector x_prev_steps = data_in( _ , 2) - data_in( _ , 4);
  NumericVector y_prev_steps = data_in( _ , 3) - data_in( _ , 5);
  NumericVector prev_headings = atan2_cpp(x_prev_steps, y_prev_steps);
  NumericVector turn_angles = headings - prev_headings;
  NumericVector t_values = data_in( _ , 6);
  
  NumericVector steplengths_0 = steplengths[t_values <= times(0) | t_values > times(3) | (t_values > times(1) & t_values <= times(2))];
  NumericVector steplengths_1 = steplengths[t_values > times(0) & t_values <= times(1)];
  NumericVector steplengths_2 = steplengths[t_values > times(2) & t_values <= times(3)];
  NumericVector turnangles_0 = turn_angles[t_values <= times(0) | t_values > times(3) | (t_values > times(1) & t_values <= times(2))];
  NumericVector turnangles_1 = turn_angles[t_values > times(0) & t_values <= times(1)];
  NumericVector turnangles_2 = turn_angles[t_values > times(2) & t_values <= times(3)];
  
  double rho_0 = mean(steplengths_0);
  double rho_1 = mean(steplengths_1) - rho_0;
  if (rho_1 < 0) {
    rho_1 = 0;
  }
  double rho_2 = mean(steplengths_2) - rho_0;
  if (rho_2 < 0) {
    rho_2 = 0;
  }
  double kappa_0 = mle_vm_kappa(turnangles_0);
  double kappa_1 = mle_vm_kappa(turnangles_1) - kappa_0;
  if (kappa_1 < 0) {
    kappa_1 = 0;
  }
  double kappa_2 = mle_vm_kappa(turnangles_2) - kappa_0;
  if (kappa_2 < 0) {
    kappa_2 = 0;
  }
  
  double NLL_rho0 = sum(steplengths_0 / rho_0 - log(1 / rho_0));
  double NLL_rho1 = sum(steplengths_1 / (rho_0 + rho_1) - log(1 / (rho_0 + rho_1)));
  double NLL_rho2 = sum(steplengths_2 / (rho_0 + rho_2) - log(1 / (rho_0 + rho_2)));
  double NLL_kappa0 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0, 0, 1)) - kappa_0 * cos(turnangles_0));
  double NLL_kappa1 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0 + kappa_1, 0, 1)) - (kappa_0 + kappa_1) * cos(turnangles_1));
  double NLL_kappa2 = sum(log(M_PI * 2 * Rf_bessel_i(kappa_0 + kappa_2, 0, 1)) - (kappa_0 + kappa_2) * cos(turnangles_2));
  
  return {NLL_rho0 + NLL_rho1 + NLL_rho2 + NLL_kappa0 + NLL_kappa1 + NLL_kappa2, rho_0, rho_1, rho_2, kappa_0, kappa_1, kappa_2};
  
}

// [[Rcpp::export]]
NumericVector fit_migration_fixed_times(NumericMatrix data_in, NumericVector times, String mod_type = "r1k1") {
  
  int n_times = times.size();
  // Rcout << "We're in \n";
  
  if ((mod_type == "r1k1" | mod_type == "rCkC") & n_times == 2) {
    // Rcout << "We're inside the if statement 1\n";
    return fit_1migration(data_in, times);
  } else if (mod_type == "r1k1") {
    // Rcout << "We're inside the if statement 2\n";
    return fit_2migrationr1k1(data_in, times);
  } else if (mod_type == "rCkC") {
    // right now we don't have capability for c > 2 but will come eventually
    // Rcout << "We're inside the if statement 3\n";
    return fit_2migrationrCkC(data_in, times);
  } else {
    // Rcout << "We're inside the if statement 4\n";
    stop("invalid mod_type argument; please select r1k1 or rCkC");
  }
  
}

// Fits the multivariate Gaussian migration model using TMB (r(t) is a step function with fixed bounds)
//
// times_list: matrix or data.frame with 2c columns; assumed to meet 0 < [i, 0] < [i, 1] < ... < [i, 2c-1] < max(data_in[,6]) for all i
// data_in: data.frame or matrix with fields [longitude/x, latitude/y, x(t-1), y(t-1), x(t-2), y(t-2), difftime/# of days from start]
// mod_type: what type of model are we running ("r1k1" if there are only two migratory "states" regardless of c)
//
// Returns a model fit
// [[Rcpp::export]]
NumericMatrix estimate_times_cpp(NumericMatrix times_matrix, NumericMatrix data_in, String mod_type = "r1k1") {

  if (mod_type != "r1k1" & mod_type != "rCkC") {
    throw std::invalid_argument("invalid mod_type argument; please select r1k1 or rCkC");
  }

  int n_unique_times = times_matrix.nrow();
  int n_breakpoints = times_matrix.ncol();
  int n_cols_all_results = n_breakpoints + 5;
  if (mod_type == "rCkC") {
    n_cols_all_results += n_breakpoints - 2;
  }

  NumericMatrix all_results(n_unique_times, n_cols_all_results);
  NumericMatrix all_results_pars(n_unique_times, n_cols_all_results - n_breakpoints);

  for (int row = 0; row < n_unique_times; row++) {
    all_results_pars(row, _ ) = fit_migration_fixed_times(data_in, times_matrix(row, _ ), mod_type);
  }
  
  for (int t = 0; t < n_breakpoints; t++) {
    all_results( _ , t) = times_matrix( _ , t);
  }
  
  for (int p = n_breakpoints; p < n_cols_all_results; p++) {
    all_results( _ , p) = all_results_pars( _ , p - n_breakpoints);
  }
  
  return all_results;

}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
