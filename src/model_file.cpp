#include "classes.hpp"

using namespace std;

void model::createKDtree () {
  
  intensivePrint ("Creating KD-trees.");
  
  trees.resize (numModelRegions);
  datKD.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    kdtree *tree  = kd_create (3);
    size_t numParams = x[r].size ();
    datKD[r].resize (numParams);

    for (size_t i=0; i<numParams; i++) {
      
      datKD[r][i] = i;
      kd_insert3 (tree, x[r][i], y[r][i], z[r][i], &datKD[r][i]);      
      
    }
    
    trees[r].push_back (tree);
  
  }  
  
}


void model::rotate () {
  
  if (angle != 0.) {
    
    intensivePrint ("Rotating.");
    rotation_matrix rotMatrix (angle, xRot, yRot, zRot);
    
    for (size_t r=0; r<numModelRegions; r++) {

      size_t numParams = x[r].size ();
      for (size_t i=0; i<numParams; i++) {
      
        rotMatrix.rotate (x[r][i], y[r][i], z[r][i], x[r][i], y[r][i], z[r][i]);
      
      }
    }    
    
  } else {
    
    intensivePrint ("No rotation found");
    
  }  
  
}


void model::findMinMaxRadius () {
  
  rMin = getRadius (x[0][0], y[0][0], z[0][0]);
  rMax = getRadius (x[0][0], y[0][0], z[0][0]);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      float rad = getRadius (x[r][i], y[r][i], z[r][i]);
      
      if (rad < rMin)
        rMin = rad;
      
      if (rad > rMax)
        rMax = rad;
      
    }
  }

  // Just set rmax to be R_EARTH.
  rMax = R_EARTH;
  
}

void model::findMinMaxCartesian () {
  
  xMin = x[0][0]; yMin = y[0][0]; zMin = z[0][0];
  xMax = x[0][0]; yMax = y[0][0]; zMax = z[0][0];
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      if (x[r][i] < xMin) {
        xMin     = x[r][i];
      }
      
      if (y[r][i] < yMin) { 
        yMin     = y[r][i];
      }
      
      if (z[r][i] < zMin) {
        zMin     = z[r][i];
      }
      
      if (x[r][i] > xMax) {
        xMax     = x[r][i];
      }
      
      if (y[r][i] > yMax) {
        yMax     = y[r][i];
      }
      
      if (z[r][i] > zMax) {
        zMax     = z[r][i];       
      }
      
    }
  }
  
}