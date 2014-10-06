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
  
  intensivePrint ("Finding convex hull.");
  
  vector<float> pMaxLoc (3, 0);
  deque<facet> facets;
  deque<facet>::reverse_iterator  sIter;
  vector<facet>::iterator fIter;    
  
  // Initialize center vector.
  vector<float> cnt = initializePoint (xCnt, yCnt, zCnt);
  
  // Initialize base line of triagle.
  vector<float> v0, v1, v2, v3;
  if (((xMax - xMin) >= (yMax - yMin)) && ((xMax - xMin) >= (zMax - zMin))) {
    v0 = initializePoint (x[xMinReg][xMinInd], y[xMinReg][xMinInd], z[xMinReg][xMinInd]);
    v1 = initializePoint (x[xMaxReg][xMaxInd], y[xMaxReg][xMaxInd], z[xMaxReg][xMaxInd]);    
  }
  
  if (((yMax - yMin) >= (xMax - xMin)) && ((yMax - yMin) >= (zMax - zMin))) {
    v0 = initializePoint (x[yMinReg][yMinInd], y[yMinReg][yMinInd], z[yMinReg][yMinInd]);
    v1 = initializePoint (x[yMaxReg][yMaxInd], y[yMaxReg][yMaxInd], z[yMaxReg][yMaxInd]);  
  }
  
  if (((zMax - zMin) >= (yMax - yMin)) && ((zMax - zMin) >= (xMax - xMin))) {
    v0 = initializePoint (x[zMinReg][zMinInd], y[zMinReg][zMinInd], z[zMinReg][zMinInd]);
    v1 = initializePoint (x[zMaxReg][zMaxInd], y[zMaxReg][zMaxInd], z[zMaxReg][zMaxInd]);    
  }
  
  // Find furtherest point from triangle base.
  vector<float> x0;
  float d2Line = 0;
  x0 = initializePoint (x[xMinReg][xMinInd], y[xMinReg][xMinInd], z[xMinReg][xMinInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[xMinReg][xMinInd], y[xMinReg][xMinInd], z[xMinReg][xMinInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  x0 = initializePoint (x[xMaxReg][xMaxInd], y[xMaxReg][xMaxInd], z[xMaxReg][xMaxInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[xMaxReg][xMaxInd], y[xMaxReg][xMaxInd], z[xMaxReg][xMaxInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  x0 = initializePoint (x[yMinReg][yMinInd], y[yMinReg][yMinInd], z[yMinReg][yMinInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[yMinReg][yMinInd], y[yMinReg][yMinInd], z[yMinReg][yMinInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  x0 = initializePoint (x[yMaxReg][yMaxInd], y[yMaxReg][yMaxInd], z[yMaxReg][yMaxInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[yMaxReg][yMaxInd], y[yMaxReg][yMaxInd], z[yMaxReg][yMaxInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  x0 = initializePoint (x[zMinReg][zMinInd], y[zMinReg][zMinInd], z[zMinReg][zMinInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[zMinReg][zMinInd], y[zMinReg][zMinInd], z[zMinReg][zMinInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  x0 = initializePoint (x[zMaxReg][zMaxInd], y[zMaxReg][zMaxInd], z[zMaxReg][zMaxInd]);
  if (distFromLine (x0, v0, v1) >= d2Line) {
    v2 = initializePoint (x[zMaxReg][zMaxInd], y[zMaxReg][zMaxInd], z[zMaxReg][zMaxInd]);
    d2Line = distFromLine (x0, v0, v1);
  }
  
  // initialize base facet, and find furthest point.
  float dstMax = 0; size_t regMax, indMax;
  facet f0 (v0, v1, v2, cnt);
  for (size_t r=0; r<numModelRegions; r++) {

    int numParams = x[r].size ();
    for (size_t i=0; i<numParams; i++) {
      
      float xLoc = x[r][i];
      float yLoc = y[r][i];
      float zLoc = z[r][i];
      
      float dst = abs (projWonV_Dist (xLoc, yLoc, zLoc, f0.n, f0.v0));
      if (dst > dstMax) {
        
        regMax = r;
        indMax = i;
        dstMax = dst;
        
      }
      
    }
  }
  
  // This is the further point from base facet.
  v3 = initializePoint (x[regMax][indMax], y[regMax][indMax], z[regMax][indMax]);
         
  ofstream myfile ("facets.txt", ios::out); 
  
  // initialize remaining facets.
  facet f1 (v1, v2, v3, cnt);
  facet f2 (v0, v2, v3, cnt);
  facet f3 (v0, v1, v3, cnt);
  
  std::vector<facet> initFacets;
  initFacets.push_back (f0); initFacets.push_back (f1); 
  initFacets.push_back (f2); initFacets.push_back (f3);
    
  // Get initial point set and push facets onto the stack.
  pointTracker pTrack (*this);
  for (fIter=initFacets.begin(); fIter!=initFacets.end(); ++fIter) {
    for (size_t r=0; r<numModelRegions; r++) {

      int numParams = x[r].size ();
      for (size_t i=0; i<numParams; i++) {

        float xLoc = x[r][i];
        float yLoc = y[r][i];
        float zLoc = z[r][i];

        float d = projWonV_Dist (xLoc, yLoc, zLoc, fIter->n, fIter->v0);

        if ((d > 0) && (pTrack.check (r, i))) {
          
          fIter->regionSet.push_back (r);
          fIter->indexSet.push_back  (i);        
          pTrack.setFalse            (r, i);  
          
        }
        
      }
    }
        
    // push the initialized facet to the stack.
    if (fIter->regionSet.size () != 0)
      dumpFacet (myfile, *fIter);      
      facets.push_back (*fIter);    
  }
  
  //// END INITIALIZATION /////
  
  vector<facet> visibleSet;
  vector<facet>::iterator vIter, iIter;
  while (not facets.empty ()) {
    
    facet fScratch = facets.back ();
    
    float dstMax=0;
    size_t regMax=0, indMax=0;
    for (size_t i=0; i<fScratch.regionSet.size (); i++) {
              
      size_t regLoc = fScratch.regionSet[i];
      size_t indLoc = fScratch.indexSet[i];
            
      float xLoc = x[regLoc][indLoc];
      float yLoc = y[regLoc][indLoc];
      float zLoc = z[regLoc][indLoc];
    
      float dst = projWonV_Dist (xLoc, yLoc, zLoc, fScratch.n, fScratch.v0);

      if (dst > dstMax) {
        
        dstMax = dst;
        regMax = regLoc;
        indMax = indLoc;
        pMaxLoc[0] = x[regLoc][indLoc];
        pMaxLoc[1] = y[regLoc][indLoc];
        pMaxLoc[2] = z[regLoc][indLoc];
        
      }        
              
    }
    
    for (sIter=facets.rbegin(); sIter!=facets.rend(); ++sIter) {
      
      // check to see if current face is a neighbour.
      bool trueNeighbour = sIter->checkNeighbour (fScratch);
      
      // if it is, check to see if the star point is visible.
      sIter->pMax.resize (3);
      if (trueNeighbour) {
        
        sIter->pMax[0] = x[regMax][indMax];
        sIter->pMax[1] = y[regMax][indMax];
        sIter->pMax[2] = z[regMax][indMax];
          
        float dst = projWonV_Dist (sIter->pMax[0], sIter->pMax[0], sIter->pMax[0], 
          sIter->n, sIter->v0);
        
        if (dst > 0) {
          visibleSet.push_back (*sIter);
          sIter->remove = true;
        }
        
      }
              
    }
    
    std::vector<edge> horizon;
    for (vIter=visibleSet.begin(); vIter!=visibleSet.end(); ++vIter) {
      for (iIter=visibleSet.begin(); iIter!=visibleSet.end(); ++iIter) {
                
        if (visibleSet.size () == 1) {
          horizon.push_back (edge (v0, v1));
          horizon.push_back (edge (v1, v2));
          horizon.push_back (edge (v2, v0));
          
        } else {
          
          cout << "BIG" << endl;
          iIter->findHorizon (*vIter);
          
        }
                      
      }            
    }

    std::vector<facet> newFacets;
    std::vector<edge>::iterator eIter;
    for (eIter=horizon.begin(); eIter!=horizon.end(); ++eIter) {
      
      facet f (eIter->v0, eIter->v1, pMaxLoc, cnt);
      newFacets.push_back (f);
    }

    pTrack.reset (*this);
    for (vIter=newFacets.begin(); vIter!=newFacets.end(); ++vIter) {
      for (size_t i=0; i<fScratch.regionSet.size (); i++) {
      
        size_t regLoc = fScratch.regionSet[i];
        size_t indLoc = fScratch.indexSet[i];

        float xLoc = x[regLoc][indLoc];
        float yLoc = y[regLoc][indLoc];
        float zLoc = z[regLoc][indLoc];
    
        float dst = projWonV_Dist (xLoc, yLoc, zLoc, fScratch.n, fScratch.v0);
      
        if ((dst > 0) && (pTrack.check (regLoc, indLoc))) {
        
          vIter->regionSet.push_back (regLoc);
          vIter->indexSet.push_back  (indLoc);        
          pTrack.setFalse            (regLoc, indLoc);  
        
        }
      
      }
    }
    
    deque<facet>::iterator sfIter;
    for (sfIter=facets.begin(); sfIter!=facets.end(); ) {
      if (sfIter->remove) {
        sfIter = facets.erase (sfIter);
      } else {
        ++sfIter;
      }
    }

    for (vIter=newFacets.begin(); vIter!=newFacets.end(); ++vIter) {
      dumpFacet (myfile, *vIter);
      cout << vIter->v0[0] << endl;
      facets.push_back (*vIter);
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
  xCnt = x[0][0]; yCnt = y[0][0]; zCnt = z[0][0];
  
  // TODO make total a global variable.
  float xSum=0, ySum=0, zSum=0;
  int total=0;
  for (size_t r=0; r<numModelRegions; r++) {
    
    int numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      if (x[r][i] < xMin) {
        xMin     = x[r][i];
        xMinReg = r;
        xMinInd = i;
      }
      
      if (y[r][i] < yMin) { 
        yMin     = y[r][i];
        yMinReg = r;
        yMinInd = i;
      }
      
      if (z[r][i] < zMin) {
        zMin     = z[r][i];
        zMinReg = r;
        zMinInd = i;
      }
      
      if (x[r][i] > xMax) {
        xMax     = x[r][i];
        xMaxReg = r;
        xMaxInd = i;
      }
      
      if (y[r][i] > yMax) {
        yMax     = y[r][i];
        yMaxReg = r;
        yMaxInd = i;
      }
      
      if (z[r][i] > zMax) {
        zMax     = z[r][i];       
        zMaxReg = r;
        zMaxInd = i;
      }
     
      xSum += x[r][i];
      ySum += y[r][i];
      zSum += z[r][i];
      
      total++;
      
    }
  }
  
  xCnt = xSum / total;
  yCnt = ySum / total;
  zCnt = zSum / total;
  
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

void model::dumpPointCloud () {
  
  if (myRank == 0) {
    ofstream myfile ("points.txt", ios::out);
  
    for (size_t r=0; r<numModelRegions; r++) {
    
      int numParams = x[r].size ();
      for (size_t i=0; i<numParams; i++) {
        myfile << x[r][i] << " " << y[r][i] << " " << z[r][i] << endl;      
      
      }
    }
  }
  
}

void model::dumpFacet (ofstream &myfile, facet &f) {
  
  if (myRank == 0) {
    
    myfile << f.v0[0] << " " << f.v0[1] << " " << f.v0[2] << endl;
    myfile << f.v1[0] << " " << f.v1[1] << " " << f.v1[2] << endl;
    myfile << f.v2[0] << " " << f.v2[1] << " " << f.v2[2] << endl;
    
  }
}


// if (myRank == 0) {
//     myfile << fMaster->v0[0] << " " << fMaster->v0[1] << " " << fMaster->v0[2] << endl;
//     myfile << fMaster->v1[0] << " " << fMaster->v1[1] << " " << fMaster->v1[2] << endl;
//     myfile << fMaster->v2[0] << " " << fMaster->v2[1] << " " << fMaster->v2[2] << endl;
// }

// bool finished = false;
// vector<facet> newFacets;
// while (not finished) {
//
//   // Reset the point tracker. All points in the index set are fair game.
//   pTrack.reset (*this);
//   // finished = true;
//
//   std::vector<facet> facets;
//   facets = allFacets;
//
//   // Loop over test facets.
//   int k=0;
//   for (fMaster=facets.rbegin(); fMaster!=facets.rend(); ++fMaster) {
//
//     k++;
//     // Reset all the region set information.
//     fMaster->reset ();
//
//     // Loop over all points in the previously found set.
//     for (size_t i=0; i<fMaster->oldRegionSet.size (); i++) {
//
//       size_t reg = fMaster->oldRegionSet[i];
//       size_t ind = fMaster->oldIndexSet[i];
//
//       float xLoc = x[reg][ind];
//       float yLoc = y[reg][ind];
//       float zLoc = z[reg][ind];
//
//       float d = projWonV_Dist (xLoc, yLoc, zLoc, fMaster->n, fMaster->v0);
//
//       if (d > 0 && pTrack.check (reg, ind)) {
//
//         fMaster->hasOutsideSet = true;
//         fMaster->regionSet.push_back (reg);
//         fMaster->indexSet.push_back  (ind);
//         pTrack.setFalse         (reg, ind);
//
//         if (d > (fMaster->dMax)) {
//
//           fMaster->dMax       = d;
//           fMaster->dMaxRegion = reg;
//           fMaster->dMaxIndex  = ind;
//
//         }
//       }
//     }
//
//     for (fNeighbour=facets.rbegin(); fNeighbour!=facets.rend(); ++fNeighbour) {
//
//       // assume we're not a neighbour, and then check that assumption.
//       bool trueNeighbour = false;
//       trueNeighbour      = fNeighbour->checkMaster (*fMaster);
//       cout << trueNeighbour
//
//       int reg = fMaster->dMaxRegion;
//       int ind = fMaster->dMaxIndex;
//
//       bool illuminated=false;
//
//       // if we are a neighbour, check if we're illuminated by the most distant point.
//       if (trueNeighbour && fMaster->hasOutsideSet) {
//
//         float d = projWonV_Dist (x[reg][ind], y[reg][ind], z[reg][ind],
//           fNeighbour->n, fNeighbour->v0);
//
//         if (d > 0)
//           illuminated = true;
//
//       }
//
//       if (illuminated) {
//
//         vector<float> vNew (3,0);
//         vNew[0] = x[reg][ind]; vNew[1] = y[reg][ind]; vNew[2] = z[reg][ind];
//         vector<edge> freeEdge = fNeighbour->findHorizon (*fMaster);
//         facet f0 (freeEdge[0].v0, freeEdge[0].v1, vNew, cnt);
//         facet f1 (freeEdge[1].v0, freeEdge[1].v1, vNew, cnt);
//
//         newFacets.push_back (f0);
//         newFacets.push_back (f1);
//
//         cout << "FOUND" << endl;
//         pointTracker newPtrack (*this);
//         for (vector<facet>::iterator newSearch=newFacets.begin(); newSearch!=newFacets.end(); newSearch++) {
//           for (size_t i=0; i<fMaster->regionSet.size (); i++) {
//
//             size_t reg = fMaster->regionSet[i];
//             size_t ind = fMaster->indexSet[i];
//
//             float xLoc = x[reg][ind];
//             float yLoc = y[reg][ind];
//             float zLoc = z[reg][ind];
//
//             float d = projWonV_Dist (xLoc, yLoc, zLoc, fMaster->n, fMaster->v0);
//
//
//             if (d > 0 && newPtrack.check (reg, ind)) {
//
//               newSearch->hasOutsideSet = true;
//               newSearch->regionSet.push_back (reg);
//               newSearch->indexSet.push_back  (ind);
//               newPtrack.setFalse         (reg, ind);
//
//             }
//           }
//         }
//
//
//         fNeighbour->remove = true;
//
//       }
//     }
//
//   if (fMaster->hasOutsideSet)
//     fMaster->remove = true;
//   }
//
//   vector<facet>::iterator fRemove;
//   for (fRemove=facets.begin(); fRemove!=facets.end(); ) {
//
//     if (fRemove->remove) {
//       fRemove = facets.erase(fRemove);
//     } else {
//       ++fRemove;
//     }
//
//   }
//   cout << "HI " << allFacets.size() <<  endl;
//
//   allFacets = facets;
//   allFacets.insert (allFacets.end (), newFacets.begin(), newFacets.end());
//   newFacets.clear ();
//
// }
