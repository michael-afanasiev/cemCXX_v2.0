#include "classes.hpp"

using namespace std;

ses3d::ses3d () {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
   
  double deg = 57.5;
  angle = deg2Rad (deg);
  xRot = 0.;
  yRot = 1.;
  zRot = 0.;
  
  readParameterFile   ();
  read                ();
  broadcast           ();
  convert2Radians     ();
  convert2Cartesian   ();
  rotate              ();
  findMinMaxRadius    ();

  if (direction == "extract")
    allocateArrays      ();
  
  if (direction == "interpolate")
    createKDtree ();
        
}

void ses3d::read () {
  
  intensivePrint ("Reading SES3D model.");
  if (myRank == 0) {
    
    readFile (col, "colatitude");
    readFile (lon, "longitude");
    readFile (rad, "radius");    
    
  }
  
  if (myRank == 0 && direction == "interpolate") {    

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
  
      readFile (vpv, "vpi");
      readFile (vph, "vpi");
      readFile (vsv, "vsv");
      readFile (vsh, "vsh");
  
    } else {
    
      error ("Symmetry system not yet implemented. See manual.");
    
    }  
  }
  
  if (myRank == 0 && taper == "true") {
    
    readFile (smooth, "smoother");
    
  }
  
}

void ses3d::convert2Radians () {
  
  vector<double>::iterator it;
  
  size_t k=0;
  for (size_t i=0; i<numModelRegions; i++) {
    
    k = 0;
    for (it=col[i].begin(); it!=col[i].end(); ++it) {
      
      double tmp = deg2Rad (*it);
      col[i][k] = tmp;
      k++;
      
    }
   
    k = 0;
    for (it=lon[i].begin(); it!=lon[i].end(); ++it) {
      
      double tmp = deg2Rad (*it);
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
  broadcast2DVector (smooth);
          
}

void ses3d::readFile (vector<vector<double>> &vec, string type) {
  
  std::string line;
  std::string fileName;
  std::vector<std::vector<double>> dummy;
  
  if (type == "colatitude")
    fileName = path + "/block_x";
  if (type == "longitude")
    fileName = path + "/block_y";
  if (type == "radius")
    fileName = path + "/block_z";    
  if (type == "rho")
    fileName = path + "/dRHO";
  if (type == "vpi")
    fileName = path + "/dVPP";
  if (type == "vpv")
    fileName = path + "/dVPP";
  if (type == "vph")
    fileName = path + "/dVPP";
  if (type == "vsv")
    fileName = path + "/dVSV";
  if (type == "vsh")
    fileName = path + "/dVSH";
  if (type == "smoother")
    fileName = path + "/smoother";
  
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
  
  // push back min/max regional radius values.
  if (type == "radius") {
    
    for (size_t i=0; i<dummy.size (); i++) {
      minRadRegion.push_back (*std::min_element (dummy[i].begin (), dummy[i].end ()));
      maxRadRegion.push_back (*std::max_element (dummy[i].begin (), dummy[i].end ()));
    }
    
  }
    
  if (type == "colatitude" || type == "longitude" || type == "radius") {
    vec.resize (numModelRegions);
    for (size_t i=0; i<dummy.size(); i++) {
  
      vec[i].reserve (dummy[i].size()-1);
      for (size_t j=0; j<dummy[i].size()-1; j++) {
        vec[i].push_back ((dummy[i][j] + dummy[i][j+1]) / 2.);
      }
    }
  } else {
    vec = dummy;
  }
  
}

void ses3d::write () {
  
  construct ();
  
  if (symSys.compare (0, 3, "tti") == 0) {

    writeFile (rho, "rho");
    writeFile (vpv, "vpv");
    writeFile (vph, "vph");
    writeFile (vsv, "vsv");
    writeFile (vsh, "vsh");    
    
  }
  
}

void ses3d::writeFile (std::vector<std::vector<double>> &vec, std::string fName) {
  
  // Function to write out a file in the ses3d region format.
  
  std::string omd = path + "/CEM";
  
  mkdir (omd.c_str (), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  std::cout << "Writing: " << fName << std::endl;
  ofstream myfile (omd + "/" + fName, ios::out);
  
  myfile << vec.size () << "\n";
  for (size_t r=0; r<vec.size (); r++) {
    myfile << vec[r].size () << "\n";
    for (size_t i=0; i<vec[r].size (); i++) {
      myfile << vec[r][i] << "\n";
    }
  }
  
  myfile.close ();
  
}

void ses3d::convert2Cartesian () {
  
  // This function converts the spherical co-ordinates of the model array into xyz coordinates, 
  // overwriting the old ones.
  
  intensivePrint ("Converting to cartesian co-ordinates.");
  
  vector<vector<double>>::iterator outer;
  vector<double>::iterator colIter, lonIter, radIter;
  
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
    
    size_t numParams = col[i].size () * lon[i].size() * rad[i].size ();
    x[i].resize (numParams);
    y[i].resize (numParams);
    z[i].resize (numParams);
    
    for (colIter=col[i].begin(); colIter!=col[i].end(); ++colIter) {      
      for (lonIter=lon[i].begin(); lonIter!=lon[i].end(); ++lonIter) {
        for (radIter=rad[i].begin(); radIter!=rad[i].end(); ++radIter) {
          
          colLonRad2xyz (x[i][k], y[i][k], z[i][k], *colIter, *lonIter, *radIter);
          k++;
          
        }
      }
    }    
  }
    
}
