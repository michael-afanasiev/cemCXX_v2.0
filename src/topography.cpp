#include "classes.hpp"

discontinuity::discontinuity () {
  
  readTopography ();
  
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
  tree = kd_create (3);
  
  intensivePrint ("Generating kd-tree (topography).");
  
#pragma omp parallel for
  for (size_t i=0; i<numTopParams; i++) {
    
    KDdat[i] = i;
    kd_insert3 (tree, colTop[i], lonTop[i], R_EARTH, &KDdat[i]);
    
  }
  
  cout << grn << "Done." << rst << endl;
  
  
}