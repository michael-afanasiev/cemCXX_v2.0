#include "classes.hpp"

mesh::mesh (exodus_file &eFile) {
  
  eFile.getXYZ (x, y, z);
  c11 = eFile.getVariable ("du1");  
  
}

void mesh::interpolate (model &mod) {
  
  size_t numNodes = x.size ();
  for (size_t i=0; i<numNodes; i++) {
    
    int region = mod.testBoundingBox (x[i], y[i], z[i]);
    
    // if (region == (-1))
    //   continue;
    
    c11[i] = region+1;
    
  }
  
}

void mesh::dump (exodus_file &eFile) {
  
  eFile.writeVariable (c11, "du1");
  
}