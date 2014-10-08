#include "classes.hpp"

using namespace std;

mesh::mesh (exodus_file &eFile) {
  
  eFile.getXYZ (x, y, z);
  c11 = eFile.getVariable ("du1");  
  connectivity = eFile.returnConnectivity ();
  nodeNumMap   = eFile.returnNodeNumMap   ();
  numNodes = eFile.numNodes;
  contains = eFile.returnContains ();
  
}

void mesh::interpolate (model &mod) {
  
  intensivePrint ("Interpolating.");
  size_t sizeConnect = connectivity.size ();
  int percent = sizeConnect / 100.;
  
  int pCount = 0;
  int pIter  = 0;
  for (size_t i=0; i<sizeConnect; i++) {
    
    if (not contains[i])
      continue;
         
    // extract node number.
    int nodeNum = connectivity[i] - 1;
    
    // find closest point [region specific].
    for (size_t r=0; r<mod.numModelRegions; r++) {

      kdres *set = kd_nearest3 (mod.trees[r], x[nodeNum], y[nodeNum], z[nodeNum]);
      void *ind  = kd_res_item_data (set);
      int point  = * (int *) ind;
      kd_res_free (set); 
      
      c11[nodeNum] = nodeNum;//mod.vsh[r][point];
      
    }
    
    pCount++;    
    if (pCount % percent == 0) {
      cout << pIter << " %\r" << flush;
      pIter++;
    }        
         
  }
  
  
}

void mesh::dump (exodus_file &eFile) {
  
  eFile.writeVariable (c11, "du1");
  
}