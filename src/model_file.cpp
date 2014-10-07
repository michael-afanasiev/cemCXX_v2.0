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

bool model::testBoundingBox (double x, double y, double z) {
  
  double colTest, lonTest, radTest;
  
  // rotate to plane.
  rotateToYZ->rotate (x, y, z, x, y, z);
  rotateToXY->rotate (x, y, z, x, y, z);
  
  // transform to spherical co-ordinate.
  xyz2ColLonRad (x, y, z, colTest, lonTest, radTest);

  // test to see if point is within rotated domain.
  if (colTest <= colMaxSearch && colTest >= colMinSearch &&
      lonTest <= lonMaxSearch && lonTest >= lonMinSearch &&
      radTest <= radMaxSearch && radTest >= radMinSearch) {
  
    return true;
    
  } else {
    
    return false;
    
  }          
  
}

void model::findBoundingBox () {
  
  double colCtr, lonCtr, radDum;
  xyz2ColLonRad (xCtr, yCtr, zCtr, colCtr, lonCtr, radDum);
  
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

  // Rotate to XY plane.
  double NINETY  = 90.;
  double angleXY = -1 * (deg2Rad (NINETY) - colCtr);
  xRot = 0; yRot = 1; zRot = 0;
  rotateToXY = new rotation_matrix (angleXY, xRot, yRot, zRot);
  for (size_t r=0; r<numModelRegions; r++) {

    size_t numParams = x[r].size ();
    for (size_t i=0; i<numParams; i++) {
    
      rotateToXY->rotate (xSearch[r][i], ySearch[r][i], zSearch[r][i], 
                          xSearch[r][i], ySearch[r][i], zSearch[r][i]);
    
    }
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
      
      double rad = getRadius (x[r][i], y[r][i], z[r][i]);
      
      if (rad < rMin)
        rMin = rad;
      
      if (rad > rMax)
        rMax = rad;
      
    }
  }

  // Just set rmax to be R_EARTH.
  rMax = R_EARTH;
  
}

void model::findMinMaxPhys () {
  
  xMin = x[0][0]; yMin = y[0][0]; zMin = z[0][0];
  xMax = x[0][0]; yMax = y[0][0]; zMax = z[0][0];
  xCtr = 0; yCtr = 0; zCtr = 0;
  
  colMin = 180.; colMax = 0.;
  lonMin = 180.; lonMax = -180.;  
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double colDum, lonDum, radDum;
      xyz2ColLonRad (x[r][i], y[r][i], z[r][i], colDum, lonDum, radDum);
      
      xCtr += x[r][i];
      yCtr += y[r][i];
      zCtr += z[r][i];
      
      if (x[r][i] < xMin)
        xMin   = x[r][i];
    
      if (y[r][i] < yMin)
        yMin   = y[r][i];
      
      if (z[r][i] < zMin) 
        zMin   = z[r][i];
      
      if (x[r][i] > xMax) 
        xMax   = x[r][i];
      
      if (y[r][i] > yMax) 
        yMax   = y[r][i];
      
      if (z[r][i] > zMax)
        zMax   = z[r][i];       
          
      if (colDum < colMin)
        colMin = colDum;
    
      if (colDum > colMax)
        colMax = colDum;
    
      if (lonDum < lonMin)
        lonMin = lonDum;
    
      if (lonDum > lonMax)
        lonMax = lonDum;

                  
    }
        
    xCtr = xCtr / (numParams);
    yCtr = yCtr / (numParams);
    zCtr = zCtr / (numParams);
            
  }
  
}

void model::findMinMaxRot () {
  
  xMinSearch = xSearch[0][0]; yMinSearch = ySearch[0][0]; zMinSearch = zSearch[0][0];
  xMaxSearch = xSearch[0][0]; yMaxSearch = ySearch[0][0]; zMaxSearch = zSearch[0][0];
  xCtrSearch = 0; yCtrSearch = 0; zCtrSearch = 0;
  
  colMinSearch = 180.; colMaxSearch = 0.;
  lonMinSearch = 180.; lonMaxSearch = -180.;  
  radMinSearch = R_EARTH; radMaxSearch = 0.;
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = xSearch[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double colDum, lonDum, radDum;
      xyz2ColLonRad (xSearch[r][i], ySearch[r][i], zSearch[r][i], colDum, lonDum, radDum);
      
      xCtrSearch += xSearch[r][i];
      yCtrSearch += ySearch[r][i];
      zCtrSearch += zSearch[r][i];
      
      if (xSearch[r][i] < xMin)
        xMinSearch   = xSearch[r][i];
    
      if (ySearch[r][i] < yMin)
        yMinSearch   = ySearch[r][i];
      
      if (zSearch[r][i] < zMin) 
        zMinSearch   = zSearch[r][i];
      
      if (xSearch[r][i] > xMax) 
        xMaxSearch   = xSearch[r][i];
      
      if (ySearch[r][i] > yMax) 
        yMaxSearch   = ySearch[r][i];
      
      if (zSearch[r][i] > zMax)
        zMaxSearch   = zSearch[r][i];       
          
      if (colDum < colMinSearch)
        colMinSearch = colDum;
    
      if (colDum > colMaxSearch)
        colMaxSearch = colDum;
    
      if (lonDum < lonMinSearch)
        lonMinSearch = lonDum;
    
      if (lonDum > lonMaxSearch)
        lonMaxSearch = lonDum;
      
      if (radDum > radMaxSearch)
        radMaxSearch = radDum;
      
      if (radDum < radMinSearch)
        radMinSearch = radDum;

                  
    }
        
    xCtrSearch = xCtrSearch / (numParams);
    yCtrSearch = yCtrSearch / (numParams);
    zCtrSearch = zCtrSearch / (numParams);
            
  }
  
}