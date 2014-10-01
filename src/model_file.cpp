#include "classes.hpp"

using namespace std;

void model::convert2Cartesian () {
  
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

void model::createKDtree () {
  
  intensivePrint ("Creating KD-trees.");
  
  trees.resize (numModelRegions);
  datKD.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    kdtree *tree  = kd_create (3);
    int numParams = x[r].size ();
    datKD[r].resize (numParams);

    for (size_t i=0; i<numParams; i++) {
      
      datKD[r][i] = i;
      kd_insert3 (tree, x[r][i], y[r][i], z[r][i], &datKD[r][i]);      
      
    }
    
    trees[r].push_back (tree);
  
  }  
  
}