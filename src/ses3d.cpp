#include "classes.hpp"

using namespace std;

ses3d::ses3d (string pathIn, string symSysIn) {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
  
  path   = pathIn;
  symSys = symSysIn;
  
  read              ();
  broadcast         ();
  convert2Cartesian ();
  createKDtree      ();
    
}

void ses3d::read () {
  
  intensivePrint ("Reading SES3D model.");
  if (myRank == 0) {
    
    readFile (col, "colatitude");
    readFile (lon, "longitude");
    readFile (rad, "radius");

    if (symSys == "tti") {
  
      readFile (rho, "rho");
      readFile (vpv, "vpv");
      readFile (vph, "vph");
      readFile (vsv, "vsv");
      readFile (vsh, "vsh");    
      
    } else if (symSys == "tti_noRho") {
  
      readFile (vpv, "vpv");
      readFile (vph, "vph");
      readFile (vsv, "vsv");
      readFile (vsh, "vsh");                
  
    } else if (symSys == "tti_isoVp") {

      readFile (rho, "rho");
      readFile (vpi, "vpi");
      readFile (vsv, "vsv");
      readFile (vsh, "vsh");                    
  
    } else if (symSys == "tti_noRho_isoVp") {
  
      readFile (vpi, "vph");
      readFile (vsv, "vsv");
      readFile (vsh, "vsh");
  
    } else {
    
      error ("Symmetry system not yet implemented. See manual.");
    
    }  
  }
}

void ses3d::broadcast () {
  
  intensivePrint ("Broadcasting arrays.");
  
  broadcastInteger (numModelRegions);
  
  broadcast2DVector (col);
  broadcast2DVector (lon);
  broadcast2DVector (rad);
  broadcast2DVector (rho);
  broadcast2DVector (vpv);
  broadcast2DVector (vph);
  broadcast2DVector (vsv);
  broadcast2DVector (vsh);
          
}

void ses3d::readFile (vector<vector<float>> &vec, string type) {
  
  std::string line;
  std::string fileName;
  std::vector<std::vector<float>> dummy;
  
  if (type == "colatitude")
    fileName = path + "/block_x";
  if (type == "longitude")
    fileName = path + "/block_y";
  if (type == "radius")
    fileName = path + "/block_z";    
  if (type == "rho")
    fileName = path + "/dRHO";
  if (type == "vpv")
    fileName = path + "/dVPP";
  if (type == "vph")
    fileName = path + "/dVPP";
  if (type == "vsv")
    fileName = path + "/dVSV";
  if (type == "vsh")
    fileName = path + "/dVSH";
  
  std::ifstream file (fileName);
    
  int l = 0;
  int r = 0;
  int intr = 0;
  int dumNumPointsInRegion = 0;
  if (file.good()) {
    
    std::cout << rst << "Reading SES3D parameter file: " << blu << fileName << rst 
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
        regionSize.push_back (dumNumPointsInRegion);
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