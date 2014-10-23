#include "classes.hpp"

using namespace std;

mesh::mesh (exodus_file &eFile) {
  
  eFileName = eFile.returnName ();
  getMinMaxDimensions();
  eFile.getXYZ (x, y, z);
  
  c11              = eFile.getVariable ("c11");  
  c12              = eFile.getVariable ("c12");  
  c13              = eFile.getVariable ("c13");  
  c14              = eFile.getVariable ("c14");  
  c15              = eFile.getVariable ("c15");  
  c16              = eFile.getVariable ("c16");  
  c22              = eFile.getVariable ("c22");  
  c23              = eFile.getVariable ("c23");  
  c24              = eFile.getVariable ("c24");  
  c25              = eFile.getVariable ("c25");  
  c26              = eFile.getVariable ("c26");  
  c33              = eFile.getVariable ("c33");  
  c34              = eFile.getVariable ("c34");  
  c35              = eFile.getVariable ("c35");  
  c36              = eFile.getVariable ("c36");  
  c44              = eFile.getVariable ("c44");  
  c45              = eFile.getVariable ("c45");  
  c46              = eFile.getVariable ("c46");  
  c55              = eFile.getVariable ("c55");  
  c56              = eFile.getVariable ("c56");  
  c66              = eFile.getVariable ("c66");  
  rho              = eFile.getVariable ("rho");    
  elv              = eFile.getVariable ("elv");
  du1              = eFile.getVariable ("du1");
  du2              = eFile.getVariable ("du2");
  
  connectivity     = eFile.returnConnectivity ();
  nodeNumMap       = eFile.returnNodeNumMap   ();
  numNodes         = eFile.numNodes;
  interpolatingSet = eFile.returnInterpolatingSet ();
  sideSetSide      = eFile.returnSideSetSide ();
  sideSetElem      = eFile.returnSideSetElem ();
  
  getSideSets ();
  
}

void mesh::interpolate (model &mod) {
  
  intensivePrint ("Interpolating.");
  size_t setSize = interpolatingSet.size ();
  int percent = (setSize / omp_get_max_threads ()) / 100.;
    
  int pCount = 0;
  int pIter  = 0;
  
#pragma omp parallel for firstprivate (pCount, pIter)
  for (size_t i=0; i<setSize; i++) {

    // extract node number.
    size_t nodeNum = interpolatingSet[i] - 1;        
   
    // use du2 as a scratch array to avoid doubly visiting points.
    if (du2[nodeNum] == 1)
      continue;
         
    // find closest point [region specific].
    for (size_t r=0; r<mod.numModelRegions; r++) {

      bool inRegion = checkInterpolatingRegion (x[nodeNum], y[nodeNum], z[nodeNum], 
                                                mod.minRadRegion[r], mod.maxRadRegion[r]);
                                                
      if (inRegion) {
      
        kdres *set = kd_nearest3 (mod.trees[r], x[nodeNum], y[nodeNum], z[nodeNum]);
        void *ind  = kd_res_item_data (set);
        int point  = * (int *) ind;
        kd_res_free (set); 
      
        elasticTensor moduli = breakdown (mod, x[nodeNum], y[nodeNum], z[nodeNum], 
                                          r, nodeNum, point);

        c11[nodeNum] = moduli.c11;
        c12[nodeNum] = moduli.c12;
        c13[nodeNum] = moduli.c13;
        c14[nodeNum] = moduli.c14;
        c15[nodeNum] = moduli.c15;
        c16[nodeNum] = moduli.c16;
        c22[nodeNum] = moduli.c22;
        c23[nodeNum] = moduli.c23;
        c24[nodeNum] = moduli.c24;
        c25[nodeNum] = moduli.c25;
        c26[nodeNum] = moduli.c26;
        c33[nodeNum] = moduli.c33;
        c34[nodeNum] = moduli.c34;
        c35[nodeNum] = moduli.c35;
        c36[nodeNum] = moduli.c36;
        c44[nodeNum] = moduli.c44;
        c45[nodeNum] = moduli.c45;
        c46[nodeNum] = moduli.c46;
        c55[nodeNum] = moduli.c55;
        c56[nodeNum] = moduli.c56;
        c66[nodeNum] = moduli.c66;
        rho[nodeNum] = moduli.rho;

        // mark that we've visited here.
        du2[nodeNum] = 1;
        
      }
      
    }
    
    if (omp_get_thread_num () == 0) {
      pCount++;
      if (pCount % percent == 0) {
        cout << pIter << " %\r" << flush;
        pIter++;
      }
    }
         
  }
  
  cout << grn << "Done." << rst << endl;
  
}

bool mesh::checkInterpolatingRegion (double &x, double &y, double &z, double minReg, 
                                     double maxReg) {
                                       
  // Check which kd-tree to extract from.                                                                              
  double radius = getRadius (x, y, z);
  double col, lon, rad;
  if        ((abs (radius - RAD_400) < TINY) && (radMax > RAD_400)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad + TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_400) < TINY) && (radMin < RAD_400)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad - TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_670) < TINY) && (radMax > RAD_670)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad + TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_670) < TINY) && (radMin < RAD_670)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad - TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_CMB) < TINY) && (radMax > RAD_CMB)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad + TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_CMB) < TINY) && (radMin < RAD_CMB)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad - TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_ICB) < TINY) && (radMax > RAD_ICB)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad + TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if ((abs (radius - RAD_ICB) < TINY) && (radMin < RAD_ICB)) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad - TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  } else if (radius > R_EARTH) {
    xyz2ColLonRad (x, y, z, col, lon, rad);
    rad = rad - TINY;
    colLonRad2xyz (x, y, z, col, lon, rad);
  }        

  if (minReg == radMin || maxReg == radMax)
    return true;
  
  if (minReg == radMax || maxReg == radMin)
    return false;
  
  if (minReg > radMin && minReg < radMax) {
    if (abs (radius - minReg) < BIGTINY && (radius - minReg) >= 0)
      return false;
  } else if (maxReg > radMin && maxReg < radMax) {
    if (abs (radius - maxReg) < BIGTINY && (radius - maxReg) >= 0)
      return true;
  }
      
  if (radius <= maxReg && radius >= minReg) {
    return true;
  } else {
    return false;
  }        
                                     
}

void mesh::createKDTree () {
  
  intensivePrint ("Creating KD-tree.");
  
  tree = kd_create (3);
  
  size_t interpolatingSetSize = interpolatingSet.size ();
  datKD.resize (interpolatingSetSize);
  
  for (size_t i=0; i<interpolatingSetSize; i++) {
    
    size_t nodeNum = interpolatingSet[i] - 1;
    datKD[i] = nodeNum;
    kd_insert3 (tree, x[nodeNum], y[nodeNum], z[nodeNum], &datKD[i]);
    
  }  
  
}

void mesh::extract (model &mod) {
  
  createKDTree ();
  
  intensivePrint ("Extracting.");
  
  // Parameters for search radius balloon.
  double searchRadius      = 1.;
  const double ONE_PERCENT = 0.01;
  const double TEN_PERCENT = 0.10;
  size_t sizeConnect       = connectivity.size ();  
  int totalLoops           = 0;
  
  // Loop over model regions.  
  for (size_t r=0; r<mod.numModelRegions; r++) {
    
    if (mod.maxRadRegion[r] <= radMin || mod.minRadRegion[r] >= radMax)
      continue;
       
    size_t numParams = mod.x[r].size ();
    
    // Initialize percentage reporting.
    int percent      = (numParams) / 100.;
    int pIter        = 0;    
            
#pragma omp parallel for firstprivate (searchRadius, r) schedule (guided)
    for (size_t i=0; i<numParams; i++) {

      double xTarget = mod.x[r][i];
      double yTarget = mod.y[r][i];
      double zTarget = mod.z[r][i];

        
      // Check if we're within the (coarse) mesh bounds.
      if (checkBoundingBox (xTarget, yTarget, zTarget)) {

        // Assume we haven't found the enclosing tet, and begin searching.
        bool found = false;
        while (not found) {
          
          // Define vector with target point.
          std::vector<double> p0 = returnVector (xTarget, yTarget, zTarget);
          
          // Get the set of nearest points to target point, with a dynamically set searchRadius.
          kdres *set = kd_nearest_range3 (tree, xTarget, yTarget, zTarget, searchRadius);

          // Initialize node and iterator numbers.
          size_t n0=0, n1=0, n2=0, n3=0;
          size_t i0=0, i1=0, i2=0, i3=0;
          
          // While we're not past the end of the set of nearest neighbours, and while we haven't
          // actually found the tet we're looking for.
          while (kd_res_end (set) == 0 && not found) {
        
            // Get the current index from the kd-tree.
            void *ind = kd_res_item_data (set);
            int point = * (int *) ind;
          
            // Loop over the entire connectivity array.
            for (size_t e=0; e<sizeConnect; e++) {
          
              // If the entry in connectivity array matches the index of the point we need 
              // (remember: all indices are node numbers - 1, due to node numbering starting at
              // one and c arrays starting at zero), figure out which element belongs to the point.
              // The else ifs here take care of the striding in the connectivity array.
              if ((connectivity[e] - 1) == point) {
      
                if        (e % numNodePerElem == 0) {
                  i0 = e+0; i1 = e+1;
                  i2 = e+2; i3 = e+3;
                } else if (e % numNodePerElem == 1) {
                  i0 = e-1; i1 = e+0;
                  i2 = e+1; i3 = e+2;
                } else if (e % numNodePerElem == 2) {
                  i0 = e-2; i1 = e-1;
                  i2 = e+0; i3 = e+1;
                } else if (e % numNodePerElem == 3) {
                  i0 = e-3; i1 = e-2; 
                  i2 = e-1; i3 = e+0;
                }         
                                   
                // Again, the actual node indices are the node numbers - 1.
                n0 = connectivity[i0]-1;
                n1 = connectivity[i1]-1;
                n2 = connectivity[i2]-1;
                n3 = connectivity[i3]-1;
                                                
                // Set up our four vectors which define the edge of the tet.
                std::vector<double> v0 = returnVector (x[n0], y[n0], z[n0]);
                std::vector<double> v1 = returnVector (x[n1], y[n1], z[n1]);
                std::vector<double> v2 = returnVector (x[n2], y[n2], z[n2]);
                std::vector<double> v3 = returnVector (x[n3], y[n3], z[n3]);
                  
                // If we're on a side set, check and project to the actual mesh if necessary.
                if (onSideSet[i0] && onSideSet[i1] && onSideSet[i3]) {
                  checkAndProject (v0, v1, v3, p0);
                } else if (onSideSet[i1] && onSideSet[i2] && onSideSet[i3]) {
                  checkAndProject (v1, v2, v3, p0);
                } else if (onSideSet[i0] && onSideSet[i2] && onSideSet[i3]) {
                  checkAndProject (v0, v2, v3, p0);
                } else if (onSideSet[i0] && onSideSet[i1] && onSideSet[i2]) {
                  checkAndProject (v0, v1, v2, p0);
                }
                  
                // Do the barycentric transform to test interpolation condition.
                double l0, l1, l2, l3;
                found = testInsideTet (v0, v1, v2, v3, p0, l0, l1, l2, l3);
                if (found == true) {
                
                  mod.c11[r][i] = interpolateTet (c11, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c12[r][i] = interpolateTet (c12, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c13[r][i] = interpolateTet (c13, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c14[r][i] = interpolateTet (c14, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c15[r][i] = interpolateTet (c15, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c16[r][i] = interpolateTet (c16, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c22[r][i] = interpolateTet (c22, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c23[r][i] = interpolateTet (c23, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c24[r][i] = interpolateTet (c24, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c25[r][i] = interpolateTet (c25, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c26[r][i] = interpolateTet (c26, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c33[r][i] = interpolateTet (c33, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c34[r][i] = interpolateTet (c34, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c35[r][i] = interpolateTet (c35, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c36[r][i] = interpolateTet (c36, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c44[r][i] = interpolateTet (c44, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c45[r][i] = interpolateTet (c45, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c46[r][i] = interpolateTet (c46, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c55[r][i] = interpolateTet (c55, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c56[r][i] = interpolateTet (c56, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.c66[r][i] = interpolateTet (c66, n0, n1, n2, n3, l0, l1, l2, l3); 
                  mod.rho[r][i] = interpolateTet (rho, n0, n1, n2, n3, l0, l1, l2, l3); 
                  
                  // Keep the search radius tight, and break out of loop.
                  searchRadius = searchRadius - searchRadius * ONE_PERCENT;
                  
                  break;               
                
                }            
              }          
            }        
      
            // Advance the result set if we haven't yet found our man.
            kd_res_next (set);
        
          }
        
          // Increase the search radius if we haven't yet found our man.
          if (not found)
            searchRadius = searchRadius + searchRadius * TEN_PERCENT;     

          if (searchRadius > 100) {
            cout << searchRadius << endl;
            double col, lon, rad;
            xyz2ColLonRad (xTarget, yTarget, zTarget, col, lon, rad);
            cout << rad2Deg(col) << ' ' << rad2Deg(lon) << ' ' << rad << endl;
          }

          // If we actually have a results set, let's free the memory needed for the next pass.
          if (kd_res_size (set) != 0)
            kd_res_free (set);
          
        }                        
      }  
      
      // Percent reporting TODO move this to a function.      
#pragma omp critical
      {
        totalLoops++;      
        if (totalLoops % percent == 0) {
          cout << pIter << " %\r" << flush;
          pIter++;
        }
      }      
    }
        
  }
  
  cout << grn << "Done." << rst << endl;
    
}

void mesh::getSideSets () {
  
  // This function creates side sets for the mesh. You can define your own side set condition with
  // if statements below.
  
  size_t connectivitySize = connectivity.size ();
  
  onSideSet.resize (connectivitySize);
  std::fill (onSideSet.begin (), onSideSet.end (), false);
  for (size_t i=0; i<connectivitySize; i++) {
    
    size_t posIter = connectivity[i] - 1;
    
    double radius = getRadius (x[posIter], y[posIter], z[posIter]);
    if ((abs (radius-radMax) < 0.5) || (abs (radius-radMin) < 0.5))
      onSideSet[i] = true;
    
  }
  
}

void mesh::checkAndProject (std::vector<double> &v0, std::vector<double> &v1, 
                            std::vector<double> &v2, std::vector<double> &p0) {
                              
  // Determines the distance from an arbitrary point which may be just outside the mesh. Then, this
  // point is modified and projected onto the border of the mesh.
                              
  if (abs (p0[0]-xMin) < CLOSE)
    p0[0] = p0[0] + CLOSE;

  if (abs (p0[1]-yMin) < CLOSE)
    p0[1] = p0[1] + CLOSE;

  if (abs (p0[2]-zMin) < CLOSE)
    p0[2] = p0[2] + CLOSE;

  // Get distance to point from origin. Make a plane from the 3 points defining a mesh edge face.
  std::vector<double> orig (3, 0);
  std::vector<double> n0 = getNormalVector (v0, v1, v2);
  double dist            = projWonV_Dist (v0, n0, orig);
  double pRadius         = getRadius (p0[0], p0[1], p0[2]);  
  
  // If this distance is given as negative, fix the normal so that it's positive.
  if (dist < 0) {
    n0[0] = n0[0] * (-1);
    n0[1] = n0[1] * (-1);
    n0[2] = n0[2] * (-1);
  }
  
  // Re-project with the new normal.
  dist = projWonV_Dist (v0, n0, orig);

  // If we're close a certain mesh edge, check to see if we're actually off the edge. If we are
  // project the point down (or up) so that it lies within the plane of the closest tet face. 
  // Tiny is necessary here to deal with small floating point errors.
  // Current settings: tiny (0.01), close (1).
  if        (pRadius > dist && (abs (pRadius-radMax) < CLOSE)) {

    double dif = abs (radMax - dist);
    p0[0] = p0[0] + (n0[0]) * (dif + TINY);
    p0[1] = p0[1] + (n0[1]) * (dif + TINY);
    p0[2] = p0[2] + (n0[2]) * (dif + TINY);

  } else if (pRadius < dist && (abs (pRadius-radMin) < CLOSE)) {

    double dif = abs (radMin - dist);
    p0[0] = p0[0] + n0[0] * (dif + TINY) * (-1);
    p0[1] = p0[1] + n0[1] * (dif + TINY) * (-1);
    p0[2] = p0[2] + n0[2] * (dif + TINY) * (-1);

  }
           
}

void mesh::interpolateTopography (discontinuity &topo) {
    
  intensivePrint ("Interpolating topography.");

  // Number of nodes in mesh chunk.
  size_t setSize = x.size ();    
#pragma omp parallel for
  for (size_t i=0; i<setSize; i++) {
      
    double col, lon, rad;            
    int nodeNum = i;

    // Get spherical co-ordinates for interpolation.
    xyz2ColLonRad (x[nodeNum], y[nodeNum], z[nodeNum], col, lon, rad);
    
    // Pull out both the topography and the crust parameters.
    kdres *setTop = kd_nearest3 (topo.elvTree,   rad2Deg (col), rad2Deg (lon), R_EARTH);        
    kdres *setCst = kd_nearest3 (topo.crustTree, rad2Deg (col), rad2Deg (lon), R_EARTH);        
    
    void *indTop  = kd_res_item_data (setTop);
    void *indCst  = kd_res_item_data (setCst);
    
    int pointTop  = * (int *) indTop;
    int pointCst  = * (int *) indCst;
    
    kd_res_free (setTop); 
    kd_res_free (setCst); 
    
    // Convert elevation to kilometers.
    elv[i] = topo.elv[pointTop] / 1000.;    

    /* The moho is defined in a weird way ( depth from sea level if in the 
    ocean, and depth from elevation if in the crust). First, convert crust 
    elevation to km, and then decided whether we're taking the sea level or
    or crustial surface as reference */    
    double referenceHeight;    
    if (elv[i] <= 0) {
      referenceHeight = R_EARTH;      
    } else {
      referenceHeight = R_EARTH + elv[i];
    }    
    double radMoho = referenceHeight - topo.dpCrust[pointCst];
    
    // If we're in crust, interpolate.
    if (rad > radMoho) {

      background_models backgroundMod;
      
      double vs1d, vp1d, rho1d;
      backgroundMod.prem_no220 (rad, vs1d, vp1d, rho1d);  
      
      double crustVsh = topo.vsCrust[pointCst];
      double crustVsv = topo.vsCrust[pointCst] - aniCorrection;    
      double crustRho = 0.2277 * crustVsh + 2.016;
      double crustVp  = 1.5399 * crustVsh + 0.840;

      double N = crustRho * crustVsh*crustVsh;
      double L = crustRho * crustVsv*crustVsv;
      double A = crustRho * crustVp*crustVp;
      double C = crustRho * crustVp*crustVp;
      double S = A - 2 * N;
      double F = A - 2 * L;
          
      c11[nodeNum] = C;
      c12[nodeNum] = F;
      c13[nodeNum] = F;
      c22[nodeNum] = A;
      c23[nodeNum] = S;
      c33[nodeNum] = A;
      c44[nodeNum] = N;
      c55[nodeNum] = L;
      c66[nodeNum] = L;
      rho[nodeNum] = crustRho;
      
      // Set the crust region flag.
      du1[nodeNum] = 1;
      
    }
    
  }
    
  cout << grn << "Done." << rst << endl;
      
}

double mesh::returnUpdate (vector<vector<double>> &vec, double &valMsh, 
                           size_t &reg, int &ind) {
                             
  // Small function that checks for the existance of a vector, and returns its value at [reg][pnt],
  // in addition to additing the value to vslMsh. If the vector does not exist, it just returns 
  // valMsh. This is the workhorse function tti model updates.
  //
  
  if (not vec.empty ()) {
    return vec[reg][ind] + valMsh;
  } else {
    return valMsh;
  }
                             
}

double mesh::SBTRKTUpdate (vector<vector<double>> &vec, double &valMsh, 
                           size_t &reg, int &ind) {
                             
  // Small function that checks for the existance of a vector, and returns its value at [reg][pnt],
  // in addition to additing the value to vslMsh. If the vector does not exist, it just returns 
  // valMsh. This is the workhorse function tti model updates.
                             
  if (not vec.empty ()) {
    return valMsh - vec[reg][ind];
  } else {
    return valMsh;
  }
                             
}

double mesh::returnUpdateAbsolute (vector<vector<double>> &vec, double &valMsh, 
                                   size_t &reg, int &ind, vector<vector<double>> &smooth) {
  
  // Small function that checks for the existance of a vector, and returns its value at [reg][pnt].
  // If it does not exist, it returns the default (in this case valmsh).
  
  if (not vec.empty ()) {
    if (not smooth.empty ()) {
      return vec[reg][ind] * smooth[reg][ind] + valMsh * (1 - smooth[reg][ind]);
    } else {
      return vec[reg][ind];
    }
  } else {
    return valMsh;
  }
  
}

double mesh::returnUpdate1d (vector<vector<double>> &vec, double &valMsh, size_t &reg, int &ind, 
                             double &val1d) {
   
  // Overload that returns the parameter added to a 1d background model.
  
  if (not vec.empty ()) {
    if (abs (vec[reg][ind]) > 0.1) {
      return val1d + vec[reg][ind];
    } else {
      return valMsh;
    }
  } else {
    return valMsh;
  }
  
}               
                        
void mesh::dump (exodus_file &eFile) {
  
  eFile.writeVariable (c11, "c11");
  eFile.writeVariable (c12, "c12");
  eFile.writeVariable (c13, "c13");
  eFile.writeVariable (c14, "c14");
  eFile.writeVariable (c15, "c15");
  eFile.writeVariable (c16, "c16");
  eFile.writeVariable (c22, "c22");
  eFile.writeVariable (c23, "c23");
  eFile.writeVariable (c24, "c24");
  eFile.writeVariable (c25, "c25");
  eFile.writeVariable (c26, "c26");
  eFile.writeVariable (c33, "c33");
  eFile.writeVariable (c34, "c34");
  eFile.writeVariable (c35, "c35");
  eFile.writeVariable (c36, "c36");
  eFile.writeVariable (c44, "c44");
  eFile.writeVariable (c45, "c45");
  eFile.writeVariable (c46, "c46");
  eFile.writeVariable (c55, "c55");
  eFile.writeVariable (c56, "c56");
  eFile.writeVariable (c66, "c66");
  eFile.writeVariable (rho, "rho");
  eFile.writeVariable (elv, "elv");
  eFile.writeVariable (du1, "du1");
  
}


elasticTensor mesh::breakdown (model &mod, double &x, double &y, double &z, 
                               size_t &region, size_t &mshInd, int &pnt) {
                        
  elasticTensor moduli;
  if (mod.symSys.compare (0, 3, "tti") == 0) {
        
    // Initialize values to mesh values.
    double vshMsh = sqrt (c44[mshInd] / rho[mshInd]);
    double vsvMsh = sqrt (c55[mshInd] / rho[mshInd]);
    double vpvMsh = sqrt (c22[mshInd] / rho[mshInd]);
    double vphMsh = sqrt (c22[mshInd] / rho[mshInd]);
    double rhoMsh = rho[mshInd];  

    double rhoNew = rhoMsh; 
    double vsvNew = vsvMsh;
    double vshNew = vshMsh;
    double vpvNew = vpvMsh;
    double vphNew = vphMsh;
  
    if (mod.interpolationType == "add_to_1d_background") {
    
      // need radius for 1d background.
      double rad = getRadius (x, y, z);
      // if (abs (rad - RAD_400) < TINY && radMin < RAD_400)
      //   rad = 5970.5;
    
      // get 1d background model.
      double vs1d, vp1d, rho1d;
      background_models backgroundMod;
      if (mod.onedBackground == "europe_model")
        backgroundMod.eumod (rad, vs1d, vp1d, rho1d);

      if (mod.onedBackground == "eumod_vpPrem_vsPremLt670")
        backgroundMod.eumod_vpPrem_vsPremLt670 (rad, vs1d, vp1d, rho1d);
    
      if (mod.onedBackground == "prem_no220")
        backgroundMod.prem_no220 (rad, vs1d, vp1d, rho1d);
            
      // get updated parameters.      
      rhoNew = returnUpdate1d (mod.rho, rhoMsh, region, pnt, rho1d);
      vsvNew = returnUpdate1d (mod.vsv, vsvMsh, region, pnt, vs1d);
      vshNew = returnUpdate1d (mod.vsh, vshMsh, region, pnt, vs1d);
      vpvNew = returnUpdate1d (mod.vpv, vpvMsh, region, pnt, vp1d);
      vphNew = returnUpdate1d (mod.vph, vphMsh, region, pnt, vp1d);
          
    }
    
    if (mod.interpolationType == "add_absolute_velocities") {
      
      // get updated parameters.      
      rhoNew = returnUpdateAbsolute (mod.rho, rhoMsh, region, pnt, mod.smooth);
      vsvNew = returnUpdateAbsolute (mod.vsv, vsvMsh, region, pnt, mod.smooth);
      vshNew = returnUpdateAbsolute (mod.vsh, vshMsh, region, pnt, mod.smooth);
      vpvNew = returnUpdateAbsolute (mod.vpv, vpvMsh, region, pnt, mod.smooth);
      vphNew = returnUpdateAbsolute (mod.vph, vphMsh, region, pnt, mod.smooth);
      
    }
    
    if (mod.interpolationType == "subtract_model") {
      
      // get updated parameters.      
      rhoNew = SBTRKTUpdate (mod.rho, rhoMsh, region, pnt);
      vsvNew = SBTRKTUpdate (mod.vsv, vsvMsh, region, pnt);
      vshNew = SBTRKTUpdate (mod.vsh, vshMsh, region, pnt);
      vpvNew = SBTRKTUpdate (mod.vpv, vpvMsh, region, pnt);
      vphNew = SBTRKTUpdate (mod.vph, vphMsh, region, pnt);
      
    }
   
    if (mod.interpolationType == "add_to_cem") {
      
      // get updated parameters.      
      rhoNew = returnUpdate (mod.rho, rhoMsh, region, pnt);
      vsvNew = returnUpdate (mod.vsv, vsvMsh, region, pnt);
      vshNew = returnUpdate (mod.vsh, vshMsh, region, pnt);
      vpvNew = returnUpdate (mod.vpv, vpvMsh, region, pnt);
      vphNew = returnUpdate (mod.vph, vphMsh, region, pnt);                        
      
    }
    
    if (mod.convert_to_1_second == "true") {
      
      attenuation atn;
      
      // need radius for 1d q model.
      double rad = getRadius (x, y, z);
      vshNew     = vshNew * atn.correctQL6 (rad);
      vsvNew     = vsvNew * atn.correctQL6 (rad);
      
    }
    
    double N = rhoNew * vshNew*vshNew;
    double L = rhoNew * vsvNew*vsvNew;
    double A = rhoNew * vphNew*vphNew;
    double C = rhoNew * vpvNew*vpvNew;
    double F = A - 2 * L;
    double S = A - 2 * N;
    
    moduli.c11 = C;
    moduli.c12 = F;
    moduli.c13 = F;
    moduli.c22 = A;
    moduli.c23 = S;
    moduli.c33 = A;
    moduli.c44 = N;
    moduli.c55 = L;
    moduli.c66 = L;
    moduli.rho = rhoNew;
        
    moduli.c14 = c14[mshInd];
    moduli.c15 = c15[mshInd];
    moduli.c16 = c16[mshInd];
    moduli.c24 = c24[mshInd];
    moduli.c25 = c25[mshInd];
    moduli.c26 = c26[mshInd];
    moduli.c34 = c34[mshInd];
    moduli.c35 = c35[mshInd];
    moduli.c36 = c36[mshInd];
    moduli.c45 = c45[mshInd];
    moduli.c46 = c46[mshInd];
    moduli.c56 = c56[mshInd];
                    
  } 
  
  return moduli;
                       
}   

void mesh::getMinMaxDimensions () {
  
  size_t found;
  found = eFileName.find ("col000-090");
  if (found != string::npos) {
    zMin = 0;
    zMax = R_EARTH;
  }
      
  found = eFileName.find ("col090-180");
  if (found != string::npos) {
    zMin = (-1) * R_EARTH;
    zMax = 0;    
  }
  
  found = eFileName.find ("lon000-090");
  if (found != string::npos) {
    xMin = 0;
    xMax = R_EARTH;
    yMin = 0;
    yMax = R_EARTH;
  }

  found = eFileName.find ("lon090-180");
  if (found != string::npos) {
    xMin = (-1) * R_EARTH;
    xMax = 0;
    yMin = 0;
    yMax = R_EARTH;
  }

  found = eFileName.find ("lon180-270");
  if (found != string::npos) {
    xMin = (-1) * R_EARTH;
    xMax = 0;
    yMin = (-1) * R_EARTH;
    yMax = 0;
  }

  found = eFileName.find ("lon270-360");
  if (found != string::npos) {
    xMin = 0;
    xMax = R_EARTH;
    yMin = (-1) * R_EARTH;
    yMax = 0;
  }

  size_t radLoc = eFileName.find_last_of ("rad");
  radMin = atof (eFileName.substr (radLoc+1, 4).c_str ());
  radMax = atof (eFileName.substr (radLoc+6, 4).c_str ());  
  
}

bool mesh::checkBoundingBox (double &x, double &y, double &z) {
  
  double rad = getRadius (x, y, z);
  
  if (x   <= xMax   && x   >= xMin &&
      y   <= yMax   && y   >= yMin &&
      z   <= zMax   && z   >= zMin &&
      rad <= radMax && rad >= radMin) {      
    return true;
  } else {
    return false;
  }
      
}
