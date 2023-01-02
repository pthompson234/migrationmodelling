template<class Type>
Type dvonmises(Type y, Type mu, Type kappa, int give_log = 0) {
  Type ret = kappa * (cos(y - mu)) - (log(2 * M_PI) + log(besselI(kappa, Type(0))));
  
  return ( give_log ? ret : exp(ret) );
}
VECTORIZE4_ttti(dvonmises)
