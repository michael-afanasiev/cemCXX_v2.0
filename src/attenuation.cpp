#include "classes.hpp"

double attenuation::QL6 (double &rad) {
  
  double Q;
  
  double x = (R_EARTH - rad) / 271.;
  
  if ((rad <= 6371) && (rad >= 6100)) {
    Q = 300. - 5370.82 * x * x + 14401.62 * x * x * x -13365.78 * 
      x * x * x * x + 4199.98 * x * x * x * x * x;
  }
  else if ( (rad <= 6100) && ( rad >= 5701) ) {
    Q = 165.;
  }
  else if ( (rad <=5701) && (rad >= 3480) ) {
    Q = 355.;
  }
  else if ( (rad <= 3480) && (rad >= 1221) ) {
    Q = 0.0; 
  }
  else {
    Q = 104.0;
  }
  
  return Q;
  
}

double attenuation::correctQL6 (double &rad) {
  
  tau_s[0] = 1.65159913;
  tau_s[1] = 13.66501919;
  tau_s[2] = 37.07774555;  
  D[0]     = 2.59203931;
  D[1]     = 2.48647256; 
  D[2]     = 0.07372733;
  
  double Q = 0;
  
  Q = QL6 ( rad );
  
  if ( Q == 0 ) {
    std::cout << "Error in setting Q value." << std::endl;
    exit (EXIT_FAILURE);
  }
  
  double freqRef = 1./30.;
  double freqMes = 1./1.;
  
  double angFreq = 2. * M_PI * freqRef;
  double tau     = 2. / ( M_PI * Q );
  
  double A = 0.0;
  double B = 0.0;
  
  for ( int i=0; i<nRelaxationMechanisms; i++ )
  {
    A = A + ( D[i] * angFreq * angFreq * tau_s[i] * tau_s[i] ) /
      ( 1 + angFreq * angFreq * tau_s[i] * tau_s[i] );
    B = B + ( D[i] * angFreq * tau_s[i] ) /
      ( 1 + angFreq * angFreq * tau_s[i] * tau_s[i] );
  }

  A = 1 + tau * A;
  B = tau * B;
  
  double factor = sqrt( ( 2 * ( A * A + B * B ) ) /
    ( ( A + sqrt ( A * A + B * B ) )) );
  
  double correctionFactor = factor - ( 1 / (M_PI * Q) * 
    log (freqRef / freqMes) );
  
  return correctionFactor;
  
}
