#include "classes.hpp"

void background_models::eumod (double &rad, double &vsv, double &vpv, double &rho) {
  
  double x = rad / R_EARTH;

  /* Stretch mantle. */
  if ( (rad <= 6371) && (rad >= 6291) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x - 0.035;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6291) && (rad >= 6191) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x - 0.035;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6191) && (rad >= 6051) ) {
    rho = 9.1790 - 5.9841   * x;
    vpv = 40.5988 - 33.5317 * x - 0.035;
    vsv = 16.8261 - 12.7527 * x - 0.065;
  }
  else if ( (rad <= 6051) && (rad >= 5971) ) {
    rho = 7.1089  - 3.8045  * x;
    vpv = 20.3926 - 12.2569 * x - 0.035;
    vsv = 8.9496  - 4.4597  * x - 0.065;
  }
  else if ( (rad <= 5971) && (rad >= 5771) ) {
    rho = 11.2494 - 8.0298  * x;
    vpv = 39.7027 - 32.6166 * x - 0.035;
    vsv = 22.3512 - 18.5856 * x - 0.12;
  }
  else if ( (rad <= 5771) && (rad >= 5701) ) {
    rho = 5.3197  - 1.4836 * x;
    vpv = 19.0957 - 9.8672 * x - 0.035;
    vsv = 9.9839  - 4.9324 * x - 0.12;
  }
  else if ( (rad <= 5701) && (rad >= 5600) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 29.2766 - 23.6026 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 22.3459 - 17.2473 * x - 2.0834 * x * x + 0.9783 * x * x * x - 0.12;
  }
  else if ( (rad <= 5600) && (rad >= 3630) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283  * x * x -3.0807   * x * x * x;
    vpv = 24.9520 - 40.4673 * x + 51.4832 * x * x - 26.6419 * x * x * x;
    vsv = 11.1671 - 13.7818 * x + 17.4575 * x * x - 9.2777  * x * x * x - 0.12;
  }
  else if ( (rad <= 3630) && (rad >= 3480) ) {
    rho = 7.9565  - 6.4761 * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 15.3891 - 5.3181 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 6.9254  + 1.4672 * x - 2.0834 * x * x + 0.9783 * x * x * x - 0.12;
  }
  else if ( (rad <= 3480) && (rad >= 1221.5) ) {
     rho = 12.5815 - 1.2638 * x - 3.6426 * x * x - 5.5281  * x * x * x;
     vpv = 11.0487 - 4.0362 * x + 4.8023 * x * x - 13.5732 * x * x * x;
     vsv = 0.0;                                                        
  }
  else if (rad <= 1221.5) {
     rho = 13.0885 - 8.8381 * x * x; 
     vpv = 11.2622 - 6.3640 * x * x; 
     vsv = 3.6678  - 4.4475 * x * x;
  }
  
}

void background_models::eumod_vpPrem_vsPremLt670 (double &rad, double &vsv, double &vpv, 
                                                  double &rho) {
  
  double x = rad / R_EARTH;

  /* Stretch mantle. */
  if ( (rad <= 6371) && (rad >= 6291) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6291) && (rad >= 6191) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6191) && (rad >= 6051) ) {
    rho = 9.1790 - 5.9841   * x;
    vpv = 40.5988 - 33.5317 * x;
    vsv = 16.8261 - 12.7527 * x - 0.065;
  }
  else if ( (rad <= 6051) && (rad >= 5971) ) {
    rho = 7.1089  - 3.8045  * x;
    vpv = 20.3926 - 12.2569 * x;
    vsv = 8.9496  - 4.4597  * x - 0.065;
  }
  else if ( (rad <= 5971) && (rad >= 5771) ) {
    rho = 11.2494 - 8.0298  * x;
    vpv = 39.7027 - 32.6166 * x;
    vsv = 22.3512 - 18.5856 * x - 0.12;
  }
  else if ( (rad <= 5771) && (rad >= 5701) ) {
    rho = 5.3197  - 1.4836 * x;
    vpv = 19.0957 - 9.8672 * x;
    vsv = 9.9839  - 4.9324 * x - 0.12;
  }
  else if ( (rad <= 5701) && (rad >= 5600) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 29.2766 - 23.6026 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 22.3459 - 17.2473 * x - 2.0834 * x * x + 0.9783 * x * x * x;
  }
  else if ( (rad <= 5600) && (rad >= 3630) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283  * x * x -3.0807   * x * x * x;
    vpv = 24.9520 - 40.4673 * x + 51.4832 * x * x - 26.6419 * x * x * x;
    vsv = 11.1671 - 13.7818 * x + 17.4575 * x * x - 9.2777  * x * x * x;
  }
  else if ( (rad <= 3630) && (rad >= 3480) ) {
    rho = 7.9565  - 6.4761 * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 15.3891 - 5.3181 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 6.9254  + 1.4672 * x - 2.0834 * x * x + 0.9783 * x * x * x;
  }
  else if ( (rad <= 3480) && (rad >= 1221.5) ) {
     rho = 12.5815 - 1.2638 * x - 3.6426 * x * x - 5.5281  * x * x * x;
     vpv = 11.0487 - 4.0362 * x + 4.8023 * x * x - 13.5732 * x * x * x;
     vsv = 0.0;                                                        
  }
  else if (rad <= 1221.5) {
     rho = 13.0885 - 8.8381 * x * x; 
     vpv = 11.2622 - 6.3640 * x * x; 
     vsv = 3.6678  - 4.4475 * x * x;
  }
    
}

void background_models::prem_no220 (double &rad, double &vsv, double &vpv, double &rho) {
  
  double x = rad / R_EARTH;

  /* Stretch mantle. */
  if ( (rad <= 6371) && (rad >= 6291) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x;
    vsv = 2.1519 + 2.3481 * x;
  }
  else if ( (rad <= 6291) && (rad >= 6191) ) {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x;
    vsv = 2.1519 + 2.3481 * x;
  }
  else if ( (rad <= 6191) && (rad >= 6051) ) {
    rho = 9.1790 - 5.9841   * x;
    vpv = 40.5988 - 33.5317 * x;
    vsv = 16.8261 - 12.7527 * x;
  }
  else if ( (rad <= 6051) && (rad >= 5971) ) {
    rho = 7.1089  - 3.8045  * x;
    vpv = 20.3926 - 12.2569 * x;
    vsv = 8.9496  - 4.4597  * x;
  }
  else if ( (rad <= 5971) && (rad >= 5771) ) {
    rho = 11.2494 - 8.0298  * x;
    vpv = 39.7027 - 32.6166 * x;
    vsv = 22.3512 - 18.5856 * x;
  }
  else if ( (rad <= 5771) && (rad >= 5701) ) {
    rho = 5.3197  - 1.4836 * x;
    vpv = 19.0957 - 9.8672 * x;
    vsv = 9.9839  - 4.9324 * x;
  }
  else if ( (rad <= 5701) && (rad >= 5600) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 29.2766 - 23.6026 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 22.3459 - 17.2473 * x - 2.0834 * x * x + 0.9783 * x * x * x;
  }
  else if ( (rad <= 5600) && (rad >= 3630) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283  * x * x -3.0807   * x * x * x;
    vpv = 24.9520 - 40.4673 * x + 51.4832 * x * x - 26.6419 * x * x * x;
    vsv = 11.1671 - 13.7818 * x + 17.4575 * x * x - 9.2777  * x * x * x;
  }
  else if ( (rad <= 3630) && (rad >= 3480) ) {
    rho = 7.9565  - 6.4761 * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 15.3891 - 5.3181 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 6.9254  + 1.4672 * x - 2.0834 * x * x + 0.9783 * x * x * x;
  }
  else if ( (rad <= 3480) && (rad >= 1221.5) ) {
     rho = 12.5815 - 1.2638 * x - 3.6426 * x * x - 5.5281  * x * x * x;
     vpv = 11.0487 - 4.0362 * x + 4.8023 * x * x - 13.5732 * x * x * x;
     vsv = 0.0;                                                        
  }
  else if (rad <= 1221.5) {
     rho = 13.0885 - 8.8381 * x * x; 
     vpv = 11.2622 - 6.3640 * x * x; 
     vsv = 3.6678  - 4.4475 * x * x;
  }
  
}

void background_models::prem (double &rad, double &vsv, double &vsh, double &vpv, 
                              double &vph, double &rho, double &eta_aniso) {
  
  double x = rad / R_EARTH;
  
  // No anisotropy by default.
  eta_aniso = 1.0;

  /* Stretch mantle. */
  if (rad > 6355.8) {
    rho = 2.6;
    vpv = 5.8;
    vsv = 3.2;
    vph = vpv;
    vsh = vsv;    
  } 
  else if ((rad <= 6355.8) && (rad >= 6346.6)) {
    rho = 2.9;
    vpv = 6.8;
    vsv = 3.9;
    vph = vpv;
    vsh = vsv;      
  } 
  else if ( (rad <= 6346.6) && (rad >= 6291) ) {
    rho       = 2.6910  + 0.6924 * x;
    vpv       = 0.8317  + 7.2180 * x;
    vph       = 3.5908  + 4.6172 * x;
    vsv       = 5.8582  - 1.4678 * x;
    vsh       = -1.0839 + 5.7176 * x;
    eta_aniso = 3.3687  - 2.4778 * x;
  } 
  else if ( (rad <= 6291) && (rad >= 6151) ) { // > 220
    rho       = 2.6910  + 0.6924 * x;
    vpv       = 0.8317  + 7.2180 * x;
    vph       = 3.5908  + 4.6172 * x;
    vsv       = 5.8582  - 1.4678 * x;
    vsh       = -1.0839 + 5.7176 * x;
    eta_aniso = 3.3687  - 2.4778 * x;
  } 
  else if ( (rad <= 6151) && (rad >= 5971) ) {
    rho = 7.1089  - 3.8045  * x;
    vpv = 20.3926 - 12.2569 * x;
    vsv = 8.9496  - 4.4597  * x;
    vph = vpv;
    vsh = vsv;
  }
  else if ( (rad <= 5971) && (rad >= 5771) ) {
    rho = 11.2494 - 8.0298  * x;
    vpv = 39.7027 - 32.6166 * x;
    vsv = 22.3512 - 18.5856 * x;
    vph = vpv;
    vsh = vsv;      
  }
  else if ( (rad <= 5771) && (rad >= 5701) ) {
    rho = 5.3197  - 1.4836 * x;
    vpv = 19.0957 - 9.8672 * x;
    vsv = 9.9839  - 4.9324 * x;
    vph = vpv;
    vsh = vsv;     
  }
  else if ( (rad <= 5701) && (rad >= 5600) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 29.2766 - 23.6026 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 22.3459 - 17.2473 * x - 2.0834 * x * x + 0.9783 * x * x * x;
    vph = vpv;
    vsh = vsv;     
  }
  else if ( (rad <= 5600) && (rad >= 3630) ) {
    rho = 7.9565  - 6.4761  * x + 5.5283  * x * x -3.0807   * x * x * x;
    vpv = 24.9520 - 40.4673 * x + 51.4832 * x * x - 26.6419 * x * x * x;
    vsv = 11.1671 - 13.7818 * x + 17.4575 * x * x - 9.2777  * x * x * x;
    vph = vpv;
    vsh = vsv;     
  }
  else if ( (rad <= 3630) && (rad >= 3480) ) {
    rho = 7.9565  - 6.4761 * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 15.3891 - 5.3181 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 6.9254  + 1.4672 * x - 2.0834 * x * x + 0.9783 * x * x * x;
    vph = vpv;
    vsh = vsv;     
  }
  else if ( (rad <= 3480) && (rad >= 1221.5) ) {
     rho = 12.5815 - 1.2638 * x - 3.6426 * x * x - 5.5281  * x * x * x;
     vpv = 11.0487 - 4.0362 * x + 4.8023 * x * x - 13.5732 * x * x * x;
     vsv = 0.0;                                                        
     vph = vpv;
     vsh = vsv;     
  }
  else if (rad <= 1221.5) {
     rho = 13.0885 - 8.8381 * x * x; 
     vpv = 11.2622 - 6.3640 * x * x; 
     vsv = 3.6678  - 4.4475 * x * x;
     vph = vpv;
     vsh = vsv;     
  }
  
}