#include "classes.hpp"

model_file::model_file (std::string testFormat, std::string testSymSys, std::string testPath) {
  
  // Test to make sure the model format makes sense.
  if (testFormat == "SES3D" || testFormat == "SPECFEM3D" || testFormat == "TERRAGRID") {
    std::cout << "Recognized model format:\t" << testFormat << std::flush << std::endl;
  } else {
    std::cout << "Model format not yet implemented. Plese choose one of:\n\tSES3D\n\tSPECFEM3D" 
      << "\n\tTERRAGRID" << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
  
  // Test to make sure the symmetry system makes sense.
  if (testSymSys == "TTI") {
    std::cout << "Recognized symmetry system:\t" << testSymSys << std::flush << std::endl;
  } else {
    std::cout << "Symmetry system not yet implemented. Plese choose one of:\n\tTTI" << std::flush 
      << std::endl;
    exit (EXIT_FAILURE);
  }
  
  modFormat = testFormat;
  symSys    = testSymSys;
  modPath   = testPath;
  
  if (modFormat == "SES3D") {
    
    readSES3D (colattitude, "colattitude");
    readSES3D (longitude, "longitude");
    readSES3D (radius, "radius");
    
    if (symSys == "TTI") {
      readSES3D (vsh, "vsh");
      readSES3D (vsh, "vsv");
      readSES3D (vpv, "vpv");
      readSES3D (vph, "vph");
      readSES3D (rho, "rho");
    }
  }
    
}