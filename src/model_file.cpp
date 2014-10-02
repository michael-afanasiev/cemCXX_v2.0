#include "classes.hpp"

using namespace std;

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

bool model::testEdge (float &x, float &y, float &z) {
  
  float dFace1 = projWonV_Dist (x, y, z, n1, p1);
  float dFace2 = projWonV_Dist (x, y, z, n2, p2);
  float dFace3 = projWonV_Dist (x, y, z, n3, p3);
  float dFace4 = projWonV_Dist (x, y, z, n4, p4);
  
  // if (MPI::COMM_WORLD.Get_rank () == 0)
  //   cout << dFace1 << " " << dFace2 << " " << dFace3 << " " << dFace4 << endl;
  
}

void model::rotate () {
  
  if (angle != 0.) {
    
    intensivePrint ("Rotating.");
    rotation_matrix rotMatrix (angle, xRot, yRot, zRot);
    
    for (size_t r=0; r<numModelRegions; r++) {

      int numParams = x[r].size ();
      for (size_t i=0; i<numParams; i++) {
      
        rotMatrix.rotate (x[r][i], y[r][i], z[r][i], x[r][i], y[r][i], z[r][i]);
      
      }
    }    
    
  } else {
    
    intensivePrint ("No rotation found");
    
  }  
  
}

void model::findConvexHull () {
  
  vector<facet>::iterator f;
  
  intensivePrint ("Finding convex hull.");
  
  float ZERO = 0.;
  
  vector<float> v0 (3,0), v1 (3,0), v2 (3,0), v3 (3,0);

  // define initial points.
  v0[0] = x[xMinI[0]][xMinI[1]];
  v0[1] = y[xMinI[0]][xMinI[1]];
  v0[2] = z[xMinI[0]][xMinI[1]];
  
  v1[0] = x[xMaxI[0]][xMaxI[1]];
  v1[1] = y[xMaxI[0]][xMaxI[1]];
  v1[2] = z[xMaxI[0]][xMaxI[1]];
  
  v2[0] = x[yMaxI[0]][yMaxI[1]];
  v2[1] = y[yMaxI[0]][yMaxI[1]];
  v2[2] = z[yMaxI[0]][yMaxI[1]];
  
  v3[0] = x[zMaxI[0]][zMaxI[1]];
  v3[1] = y[zMaxI[0]][zMaxI[1]];
  v3[2] = z[zMaxI[0]][zMaxI[1]];
  
  // initialize facets.
  facet f0 (v0, v1, v2);
  facet f1 (v0, v1, v3);
  facet f2 (v0, v2, v3);
  
  vector<facet> facets;
  facets.push_back (f0); facets.push_back (f1); facets.push_back (f2);
  
  for (f=facets.begin(); f!=facets.end(); ++f) {
    for (size_t r=0; r<numModelRegions; r++) {

      int numParams = x[r].size ();
      for (size_t i=0; i<numParams; i++) {

        float xLoc = x[r][i];
        float yLoc = y[r][i];
        float zLoc = z[r][i];

        float d = projWonV_Dist (xLoc, yLoc, zLoc, f->n, f->v0);

        if (d > 0) {
          f->region.push_back (r);
          f->index.push_back  (i);
        }

      }
    }
    
  }
  
  for (f=facets.begin(); f!=facets.end(); ++f) {
  
    for (size_t i=0; i<f->region.size (); i++) {

      size_t reg = f -> region[i];
      size_t ind = f -> index[i];
      
      float xLoc = x[reg][ind];
      float yLoc = y[reg][ind];
      float zLoc = z[reg][ind];

      float d = projWonV_Dist (xLoc, yLoc, zLoc, f->n, f->v0);
      
    }
      
  }
      
}

void model::findMinMaxRadius () {
  
  rMin = getRadius (x[0][0], y[0][0], z[0][0]);
  rMax = getRadius (x[0][0], y[0][0], z[0][0]);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    int numParams = x[r].size();
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
    
    int numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      if (x[r][i] < xMin) {
        xMin     = x[r][i];
        xMinI[0] = r;
        xMinI[1] = i;
      }
      
      if (y[r][i] < yMin) { 
        yMin     = y[r][i];
        yMinI[0] = r;
        yMinI[1] = i;
      }
      
      if (z[r][i] < zMin) {
        zMin     = z[r][i];
        zMinI[0] = r;
        zMinI[1] = i;
      }
      
      if (x[r][i] > xMax) {
        xMax     = x[r][i];
        xMaxI[0] = r;
        xMaxI[1] = i;
      }
      
      if (y[r][i] > yMax) {
        yMax     = y[r][i];
        yMaxI[0] = r;
        yMaxI[1] = i;
      }
      
      if (z[r][i] > zMax) {
        zMax     = z[r][i];       
        yMaxI[0] = r;
        yMaxI[1] = i;
      }
      
    }
  }
  
}

void model::findEdgePlanes () {
  
  intensivePrint ("Finding edges.");
  
  // Re-usable vectors for edges of chunk (for face plane).  
  std::vector<float> A, B, C;
  A.resize (3);
  B.resize (3);
  C.resize (3);

  // Normal face plane vectors.
  n1.resize (3);
  n2.resize (3);
  n3.resize (3);
  n4.resize (3);
  
  // Save edges of chunks for use to define plane later.
  p1.resize (3);
  p2.resize (3);
  p3.resize (3);
  p4.resize (3);
  
  // Face 1.
  A[0] = xMax; A[1] = yMin; A[2] = zMax;
  B[0] = xMin; B[1] = yMin; B[2] = zMax;
  C[0] = xMax; C[1] = yMin; C[2] = zMin;
  n1 = getNormalVector (A, B, C);
  p1 = A;

  // Face 2.
  A[0] = xMax; A[1] = yMax; A[2] = zMax;
  B[0] = xMax; B[1] = yMin; B[2] = zMax;
  C[0] = xMax; C[1] = yMax; C[2] = zMin;
  n2 = getNormalVector (A, B, C);
  p2 = A;

  // Face 3.
  A[0] = xMin; A[1] = yMax; A[2] = zMax;
  B[0] = xMax; B[1] = yMax; B[2] = zMax;
  C[0] = xMin; C[1] = yMax; C[2] = zMin;  
  n3 = getNormalVector (A, B, C);
  p3 = A;

  // Face 4.
  A[0] = xMin; A[1] = yMin; A[2] = zMax;
  B[0] = xMin; B[1] = yMax; B[2] = zMax;
  C[0] = xMin; C[1] = yMin; C[2] = zMin;
  n4 = getNormalVector (A, B, C);
  p4 = A;      
  
}