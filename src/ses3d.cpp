#include "classes.hpp"

void model_file::readSES3D (std::vector<std::vector<float>> &vec, std::string type) {
  
  std::string line;
  std::string fileName;
  std::vector<std::vector<float>> dummy;
  
  if (type == "colattitude")
    fileName = modPath + "/block_x";
  if (type == "longitude")
    fileName = modPath + "/block_y";
  if (type == "radius")
    fileName = modPath + "/block_z";    
  if (type == "rho")
    fileName = modPath + "/dRHO";
  if (type == "vpv")
    fileName = modPath + "/dVPP";
  if (type == "vph")
    fileName = modPath + "/dVPP";
  if (type == "vsv")
    fileName = modPath + "/dVSV";
  if (type == "vsh")
    fileName = modPath + "/dVSH";
  
  std::ifstream file (fileName);
    
  int l = 0;
  int r = 0;
  int intr = 0;
  int dumNumPointsInRegion = 0;
  int numModelRegions   = 0;
  if (file.good()) {
    
    std::cout << "Reading SES3D parameter file: " << blu << fileName << rst 
      << std::flush << std::endl;
    
    while (getline (file, line)) {
      
      // The first line is special. Get the number of ses3d regions.
      if (l == 0) {
        numModelRegions = atoi (line.c_str());
        dummy.resize (numModelRegions);
        l++;
        continue;
      }
            
      if (l == 1 || intr == dumNumPointsInRegion) {
        dumNumPointsInRegion = atoi (line.c_str());   
        dummy[r].reserve (dumNumPointsInRegion);
        intr = 0;
        l++;
        r++;
        continue;
      }
      
      if (l > 1) {
        dummy[r-1].push_back (stof (line));
        intr++;
        l++;
      }  
    }
    
  } else {
    
    std::cout << red << "Problem reading file: " << fileName << std::flush << std::endl;
    exit (EXIT_FAILURE);
    
  }
  
  int coordinateFixer;
  if (type == "colatitude" || type == "longitude" || type == "radius") {
    coordinateFixer = 1;
  } else {
    coordinateFixer = 0;
  }
    
  vec.resize (numModelRegions);
  for (size_t i=0; i<dummy.size(); i++) {
  
    vec[i].reserve (dummy[i].size()-coordinateFixer);
    for (size_t j=0; j<dummy[i].size()-coordinateFixer; j++) {
      vec[i].push_back ((dummy[i][j] + dummy[i][j+1]) / 2.);
    }
  }
  
}