#include "classes.hpp"

using namespace std;

void model::createKDtree () {
  
  intensivePrint ("Creating KD-trees.");
  
  trees.reserve (numModelRegions);
  datKD.resize  (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    kdtree *tree  = kd_create (3);
    size_t numParams = x[r].size ();
    datKD[r].resize (numParams);

    for (size_t i=0; i<numParams; i++) {
      
      datKD[r][i] = i;
      kd_insert3 (tree, x[r][i], y[r][i], z[r][i], &datKD[r][i]);      
      
    }
    
    trees.push_back (tree);
  
  }  
  
}

int model::testBoundingBox (double x, double y, double z) {
  
  int region = -1;
  double colTest, lonTest, radTest;
  
  // rotate to plane.
  rotateToYZ->rotate (x, y, z, x, y, z);
  // rotateToXY->rotate (x, y, z, x, y, z);
  
  // transform to spherical co-ordinate.
  xyz2ColLonRad (x, y, z, colTest, lonTest, radTest);

  // test to see if point is within rotated domain.
  for (size_t r=0; r<numModelRegions; r++) {
    if (colTest <= colMaxSearch[r] && colTest >= colMinSearch[r] &&
        lonTest <= lonMaxSearch[r] && lonTest >= lonMinSearch[r] &&
        radTest <= radMaxSearch[r] && radTest >= radMinSearch[r])  
          
      region = r;    
    
  }
  
  return region;
  
}

void model::findBoundingBox () {
  
  // find the vector pointing to the center of the box.
  double colCtr, lonCtr, radDum;
  xyz2ColLonRad (xCtr, yCtr, zCtr, colCtr, lonCtr, radDum);
  
  // initialize arrays to physical domain.
  xSearch = x;
  ySearch = y;
  zSearch = z;
  
  // rotate to YZ plane.
  double angleYZ = lonCtr;
  double xRot, yRot, zRot;
  xRot = 0; yRot = 0; zRot = 1;
  rotateToYZ = new rotation_matrix (angleYZ, xRot, yRot, zRot);
  
  for (size_t r=0; r<numModelRegions; r++) {

    size_t numParams = x[r].size ();
    for (size_t i=0; i<numParams; i++) {
    
      rotateToYZ->rotate (xSearch[r][i], ySearch[r][i], zSearch[r][i], 
                          xSearch[r][i], ySearch[r][i], zSearch[r][i]);
    
    }
  }   

  // // Rotate to XY plane.
  // double NINETY  = 90.;
  // double angleXY = -1 * (deg2Rad (NINETY) - colCtr);
  // xRot = 0; yRot = 1; zRot = 0;
  // rotateToXY = new rotation_matrix (angleXY, xRot, yRot, zRot);
  // for (size_t r=0; r<numModelRegions; r++) {
  //
  //   size_t numParams = x[r].size ();
  //   for (size_t i=0; i<numParams; i++) {
  //
  //     rotateToXY->rotate (xSearch[r][i], ySearch[r][i], zSearch[r][i],
  //                         xSearch[r][i], ySearch[r][i], zSearch[r][i]);
  //
  //   }
  // }
    
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
  
  rMin.resize (numModelRegions);
  rMax.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    rMin[r] = getRadius (x[r][0], y[r][0], z[r][0]);
    rMax[r] = getRadius (x[r][0], y[r][0], z[r][0]);
  }
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double rad = getRadius (x[r][i], y[r][i], z[r][i]);
      
      if (rad < rMin[r])
        rMin[r] = rad;
      
      if (rad > rMax[r])
        rMax[r] = rad;
      
    }
  }
  
}

void model::findMinMaxPhys () {
  
  xMin.resize (numModelRegions); 
  yMin.resize (numModelRegions);
  zMin.resize (numModelRegions);
  xMax.resize (numModelRegions); 
  yMax.resize (numModelRegions); 
  zMax.resize (numModelRegions);
  colMin.resize (numModelRegions); 
  colMax.resize (numModelRegions);
  lonMin.resize (numModelRegions); 
  lonMax.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    xMin[r] = x[r][0]; yMin[r] = y[r][0]; zMin[r] = z[r][0];
    xMax[r] = x[r][0]; yMax[r] = y[r][0]; zMax[r] = z[r][0];
    colMin[r] = 180.;  colMax[r] = 0.;
    lonMin[r] = 180.;  lonMax[r] = -180.;  
  }
  
  xCtr = 0; yCtr = 0; zCtr = 0;    
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double colDum, lonDum, radDum;
      xyz2ColLonRad (x[r][i], y[r][i], z[r][i], colDum, lonDum, radDum);
      
      xCtr += x[r][i];
      yCtr += y[r][i];
      zCtr += z[r][i];
      
      if (x[r][i] < xMin[r])
        xMin[r]   = x[r][i];
    
      if (y[r][i] < yMin[r])
        yMin[r]   = y[r][i];
      
      if (z[r][i] < zMin[r]) 
        zMin[r]   = z[r][i];
      
      if (x[r][i] > xMax[r]) 
        xMax[r]   = x[r][i];
      
      if (y[r][i] > yMax[r]) 
        yMax[r]   = y[r][i];
      
      if (z[r][i] > zMax[r])
        zMax[r]   = z[r][i];       
          
      if (colDum < colMin[r])
        colMin[r] = colDum;
    
      if (colDum > colMax[r])
        colMax[r] = colDum;
    
      if (lonDum < lonMin[r])
        lonMin[r] = lonDum;
    
      if (lonDum > lonMax[r])
        lonMax[r] = lonDum;

                  
    }
        
    xCtr = xCtr / (numParams);
    yCtr = yCtr / (numParams);
    zCtr = zCtr / (numParams);
            
  }
  
}

void model::findMinMaxRot () {
  
  xMinSearch.resize (numModelRegions); 
  yMinSearch.resize (numModelRegions); 
  zMinSearch.resize (numModelRegions);
  xMaxSearch.resize (numModelRegions); 
  yMaxSearch.resize (numModelRegions); 
  zMaxSearch.resize (numModelRegions);
  colMinSearch.resize (numModelRegions); 
  colMaxSearch.resize (numModelRegions);
  lonMinSearch.resize (numModelRegions); 
  lonMaxSearch.resize (numModelRegions);
  radMaxSearch.resize (numModelRegions);
  radMinSearch.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    xMinSearch[r] = x[r][0]; yMinSearch[r] = y[r][0]; zMinSearch[r] = z[r][0];
    xMaxSearch[r] = x[r][0]; yMaxSearch[r] = y[r][0]; zMaxSearch[r] = z[r][0];
    colMinSearch[r] = 180.;  colMaxSearch[r] = 0.;
    lonMinSearch[r] = 180.;  lonMaxSearch[r] = -180.;  
  }
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = xSearch[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double colDum, lonDum, radDum;
      xyz2ColLonRad (xSearch[r][i], ySearch[r][i], zSearch[r][i], colDum, lonDum, radDum);
      
      xCtrSearch += xSearch[r][i];
      yCtrSearch += ySearch[r][i];
      zCtrSearch += zSearch[r][i];
      
      if (xSearch[r][i] < xMinSearch[r])
        xMinSearch[r]   = xSearch[r][i];
    
      if (ySearch[r][i] < yMinSearch[r])
        yMinSearch[r]   = ySearch[r][i];
      
      if (zSearch[r][i] < zMinSearch[r]) 
        zMinSearch[r]   = zSearch[r][i];
      
      if (xSearch[r][i] > xMaxSearch[r]) 
        xMaxSearch[r]   = xSearch[r][i];
      
      if (ySearch[r][i] > yMaxSearch[r]) 
        yMaxSearch[r]   = ySearch[r][i];
      
      if (zSearch[r][i] > zMaxSearch[r])
        zMaxSearch[r]   = zSearch[r][i];       
          
      if (colDum < colMinSearch[r])
        colMinSearch[r] = colDum;
    
      if (colDum > colMaxSearch[r])
        colMaxSearch[r] = colDum;
    
      if (lonDum < lonMinSearch[r])
        lonMinSearch[r] = lonDum;
    
      if (lonDum > lonMaxSearch[r])
        lonMaxSearch[r] = lonDum;
      
      if (radDum > radMaxSearch[r])
        radMaxSearch[r] = radDum;
      
      if (radDum < radMinSearch[r])
        radMinSearch[r] = radDum;

                  
    }
        
    xCtrSearch = xCtrSearch / (numParams);
    yCtrSearch = yCtrSearch / (numParams);
    zCtrSearch = zCtrSearch / (numParams);
            
  }
  
}