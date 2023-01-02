#include "TMB.hpp"
#include "distributions_R_new.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( steplengths ); // step lengths
  DATA_VECTOR( turn_angles ); // turning angles (from -pi to pi)
  DATA_VECTOR( time_indices ); // number of days from the beginning of the time series for each step
  DATA_SCALAR( t1 ); // start of first migration
  DATA_SCALAR( t2 ); // end of first migration
  
  // Parameters
  PARAMETER( r0 );
  PARAMETER( r1 );
  PARAMETER( k0 );
  PARAMETER( k1 );
  
  Type val = 0.0; // negative log-likelihood
  int n_pts = time_indices.size(); // number of steps
  
  vector<Type> ll_step(n_pts);
  vector<Type> ll_turn(n_pts);
  
  for (int i = 0; i < n_pts; i++) {
    
    if (time_indices(i) > t1 && time_indices(i) < t2) {
      ll_step(i) = dexp(steplengths(i), 1.0 / (r0 + r1), true);
      ll_turn(i) = dvonmises(turn_angles(i), Type(0), k0 + k1, true);
    } else {
      ll_step(i) = dexp(steplengths(i), 1.0 / r0, true);
      ll_turn(i) = dvonmises(turn_angles(i), Type(0), k0, true);
    }

  }
  
  val = -ll_step.sum() - ll_turn.sum();
  
  // Reporting
  REPORT(val);
  REPORT( r0 );
  REPORT( r1 );
  REPORT( k0 );
  REPORT( k1 );
  
  ADREPORT( r0 );
  ADREPORT( r1 );
  ADREPORT( k0 );
  ADREPORT( k1 );
  
  return val;
}
