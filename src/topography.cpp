#include "classes.hpp"

discontinuity::discontinuity () {
  
  readTopography  ();
  readCrust       ();
  
  crustTree = createKDTree (xCrust, yCrust, crustDat);
  
}

void discontinuity::readCrust () {
  
  readCrustFile (xCrust,  "./mod/data/crust_x_smooth");
  readCrustFile (yCrust,  "./mod/data/crust_y_smooth");
  readCrustFile (vsCrust, "./mod/data/crust_vs_smooth");
  readCrustFile (dpCrust, "./mod/data/crust_dep_smooth");
  
}

void discontinuity::readCrustFile (std::vector<double> &vec, std::string fName) {

  std::ifstream myFile (fName);
  
  std::cout << "Reading crust files." << std::endl;
  if (myFile.good ()) {
    
    std::string line;
    while (std::getline (myFile, line)) {
      
      vec.push_back (stod (line));
      
    }
    
  } else {
    
    error ("Crust file " + fName + " does not exist.");
    
  }
  
}

kdtree *discontinuity::createKDTree (std::vector<double> &x, std::vector<double> &y,
                                          std::vector<int> &kdVec) {
  
  intensivePrint ("Generating kdtree (crust).");
  
  size_t size = x.size () * y.size ();
  
  kdtree *tree;
  tree = kd_create (3);
  kdVec.resize (size);
  
  int k = 0;
  for (size_t i=0; i<x.size (); i++) {
    for (size_t j=0; j<y.size (); j++) {

      if (y[j] > 180.)
        y[j] = y[j] - 360.;
          
      kdVec[k] = k;
      kd_insert3 (tree, x[i], y[j], R_EARTH, &kdVec[k]);    
    
      k++;
      
    }    
  }
  
  return tree;
  
}




void discontinuity::readTopography () {
  
  std::ifstream myFile ("./mod/data/10MinuteTopoGrid.txt");
  
  std::cout << "Reading topography." << std::endl;
  
  std::vector<double> lonTop, colTop;
  if (myFile.good ()) {
    
    std::string line;
    while (std::getline (myFile, line)) {
      
      std::stringstream linestream (line);
      std::string value;
      
      int i = 0;
      while (std::getline (linestream, value, ',')) {
        
        if (i == 0) 
          lonTop.push_back (stod (value));
        if (i == 1)
          colTop.push_back (90. - stod (value));
        if (i == 2)
          elv.push_back (stod (value));
        
        i++;
          
      }
    }    
    
  } else {
    
    error ("I can't read the topography file (should be mod/data/10MinuteTopoGrid.txt).");
    
  }
  
  size_t numTopParams = lonTop.size ();
  KDdat.resize (numTopParams);
  elvTree = kd_create (3);
  
  intensivePrint ("Generating kd-tree (topography).");
  
#pragma omp parallel for
  for (size_t i=0; i<numTopParams; i++) {
    
    KDdat[i] = i;
    kd_insert3 (elvTree, colTop[i], lonTop[i], R_EARTH, &KDdat[i]);
    
  }
  
  cout << grn << "Done." << rst << endl;
    
}
