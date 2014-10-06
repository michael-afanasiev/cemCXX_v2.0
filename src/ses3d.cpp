#include "classes.hpp"

using namespace std;

ses3d::ses3d (string pathIn, string symSysIn) {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
  
  path   = pathIn;
  symSys = symSysIn;
  
  float deg=0;
  angle = deg2Rad (deg);
  xRot = 0.;
  yRot = 1.;
  zRot = 0.;
  
  read              ();
  broadcast         ();
  convert2Radians   ();
  convert2Cartesian ();
  rotate            ();
  getTotalParameters ();
  findMinMaxCartesian        ();
  findMinMaxRadius ();
  dumpPointCloud ();
  findConvexHull ();
  findEdgePlanes    ();
  // createKDtree      ();
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    int numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      testEdge (x[r][i], y[r][i], z[r][i]);
      
    }
  }
    
}

int ses3d::getTotalParameters () {
  
  int numParams = 0;
  for (size_t i=0; i<numModelRegions; i++) {
    
    numParams += x[i].size();
    
  }
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

void ses3d::convert2Radians () {
  
  vector<float>::iterator it;
  
  size_t k=0;
  for (size_t i=0; i<numModelRegions; i++) {
    
    k = 0;
    for (it=col[i].begin(); it!=col[i].end(); ++it) {
      
      float tmp = deg2Rad (*it);
      col[i][k] = tmp;
      k++;
      
    }
   
    k = 0;
    for (it=lon[i].begin(); it!=lon[i].end(); ++it) {
      
      float tmp = deg2Rad (*it);
      lon[i][k] = tmp;
      k++;
      
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

void ses3d::convert2Cartesian () {
  
  // This function converts the spherical co-ordinates of the model array into xyz coordinates, 
  // overwriting the old ones.
  
  intensivePrint ("Converting to cartesian co-ordinates.");
  
  vector<vector<float>>::iterator outer;
  vector<float>::iterator colIter, lonIter, radIter;
  
  if (col.empty ())
    error ("No spherical co-ordinate arrays stored. Are you sure you read them in?");
  
  if (not x.empty ()) {    
    x.clear ();
    y.clear ();
    z.clear ();    
  }
  
  x.resize (numModelRegions);
  y.resize (numModelRegions);
  z.resize (numModelRegions);
  
  for (size_t i=0; i<numModelRegions; i++) {
    
    size_t k=0;
    
    int numParams = col[i].size () * lon[i].size() * rad[i].size ();
    x[i].resize (numParams);
    y[i].resize (numParams);
    z[i].resize (numParams);
    
    for (colIter=col[i].begin(); colIter!=col[i].end(); ++colIter) {
      for (lonIter=lon[i].begin(); lonIter!=lon[i].end(); ++lonIter) {
        for (radIter=rad[i].begin(); radIter!=rad[i].end(); ++radIter) {
          
          x[i][k] = *radIter * cos (*lonIter) * sin (*colIter);
          y[i][k] = *radIter * sin (*lonIter) * sin (*colIter);
          z[i][k] = *radIter * cos (*colIter);
          k++;
          
        }
      }
    }    
  }
    
}