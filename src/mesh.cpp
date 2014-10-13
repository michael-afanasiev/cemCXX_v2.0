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
  int percent = setSize / 100.;
    
  int pCount = 0;
  int pIter  = 0;
  for (size_t i=0; i<setSize; i++) {
         
    // extract node number.
    size_t nodeNum = interpolatingSet[i] - 1;
    
    // find closest point [region specific].
    for (size_t r=0; r<mod.numModelRegions; r++) {

      kdres *set = kd_nearest3 (mod.trees[r], x[nodeNum], y[nodeNum], z[nodeNum]);
      void *ind  = kd_res_item_data (set);
      int point  = * (int *) ind;
      kd_res_free (set); 
      
      elasticTensor moduli = breakdown (mod, x[nodeNum], y[nodeNum], z[nodeNum], r, nodeNum, point);
      c66[nodeNum] = moduli.c66;
      
    }
    
    pCount++;    
    if (pCount % percent == 0) {
      cout << pIter << " %\r" << flush;
      pIter++;
    }        
         
  }
  
}

void mesh::createKDTree () {
  
  intensivePrint ("Creating KD-tree.");
  datKD.resize (numNodes);
  
  tree = kd_create (3);
  
  for (size_t i=0; i<numNodes; i++) {
    
    datKD[i] = i;
    kd_insert3 (tree, x[i], y[i], z[i], &datKD[i]);
    
  }  
  
}

void mesh::extract (model &mod) {
  
  createKDTree ();
  
  intensivePrint ("Extracting.");
  
  double searchRadius = 1.;
  const double ONE_PERCENT=0.01;
  const double TEN_PERCENT=0.10;
  size_t sizeConnect = connectivity.size ();
  
  
  for (size_t r=0; r<mod.numModelRegions; r++) {
    
    size_t numParams = mod.x[r].size ();
    int percent = (numParams) / 100.;
    int pCount = 0;
    int pIter  = 0;    
        
#pragma omp parallel for firstprivate (searchRadius, pCount, pIter)
    for (size_t i=0; i<numParams; i++) {

      double xTarget = mod.x[r][i];
      double yTarget = mod.y[r][i];
      double zTarget = mod.z[r][i];
      
      if (checkBoundingBox (xTarget, yTarget, zTarget)) {
        
        bool found = false;
        while (not found) {
          
          kdres *set;             
          std::vector<double> p0 = returnVector (xTarget, yTarget, zTarget);
          
          set = kd_nearest_range3 (tree, xTarget, yTarget, zTarget, searchRadius);

          size_t n0=0, n1=0, n2=0, n3=0;
          size_t i0=0, i1=0, i2=0, i3=0;
          while (kd_res_end (set) == 0 && not found) {
        
            void *ind = kd_res_item_data (set);
            int point = * (int *) ind;
          
            for (size_t e=0; e<sizeConnect; e++) {
          
              if ((connectivity[e] - 1) == point) {
      
                if        (e % numNodePerElem == 0) {
                  i0 = e+0;
                  i1 = e+1;
                  i2 = e+2;
                  i3 = e+3;
                } else if (e % numNodePerElem == 1) {
                  i0 = e-1;
                  i1 = e+0;
                  i2 = e+1;
                  i3 = e+2;
                } else if (e % numNodePerElem == 2) {
                  i0 = e-2;
                  i1 = e-1;
                  i2 = e+0;
                  i3 = e+1;
                } else if (e % numNodePerElem == 3) {
                  i0 = e-3;
                  i1 = e-2;
                  i2 = e-1;
                  i3 = e+0;
                }         
                                   
                n0 = connectivity[i0]-1;
                n1 = connectivity[i1]-1;
                n2 = connectivity[i2]-1;
                n3 = connectivity[i3]-1;
                                                
                std::vector<double> v0 = returnVector (x[n0], y[n0], z[n0]);
                std::vector<double> v1 = returnVector (x[n1], y[n1], z[n1]);
                std::vector<double> v2 = returnVector (x[n2], y[n2], z[n2]);
                std::vector<double> v3 = returnVector (x[n3], y[n3], z[n3]);
                  
                if (onSideSet[i0] && onSideSet[i1] && onSideSet[i3]) {
                  checkAndProject (v0, v1, v3, p0);
                } else if (onSideSet[i1] && onSideSet[i2] && onSideSet[i3]) {
                  checkAndProject (v1, v2, v3, p0);
                } else if (onSideSet[i0] && onSideSet[i2] && onSideSet[i3]) {
                  checkAndProject (v0, v2, v3, p0);
                } else if (onSideSet[i0] && onSideSet[i1] && onSideSet[i2]) {
                  checkAndProject (v0, v1, v2, p0);
                }
                  
                double l0, l1, l2, l3;
                found = testInsideTet (v0, v1, v2, v3, p0, l0, l1, l2, l3);
                // cout << l0 << ' ' << l1 << ' ' << l2 << ' ' << l3 << endl;
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
                  mod.rho[r][i] = getRadius (p0[0], p0[1], p0[2]);//interpolateTet (rho, n0, n1, n2, n3, l0, l1, l2, l3); 
                  
                  searchRadius = searchRadius - searchRadius * ONE_PERCENT;
                  break;               
                
                }            
              }          
            }        
      
            kd_res_next (set);
        
          }
        
          if (not found)
            searchRadius = searchRadius + searchRadius * TEN_PERCENT;              

          if (kd_res_size (set) != 0)
            kd_res_free (set);
          
          // cout << searchRadius << ' ' << kd_res_size(set) << endl;
          // cout << getRadius (xTarget, yTarget, zTarget) << endl;
          
        }                        
      }
  
      pCount++;
      if (pCount % percent == 0) {
        cout << pIter << " %\r" << flush;
        pIter++;
      }
      // cout << numParams - i << endl;
      
    }
        
  }
    
}

void mesh::getSideSets () {
  
  size_t connectivitySize = connectivity.size ();
  
  onSideSet.resize (connectivitySize);
  std::fill (onSideSet.begin (), onSideSet.end (), false);
  for (size_t i=0; i<connectivitySize; i++) {
    
    size_t posIter = connectivity[i] - 1;
    
    double radius = getRadius (x[posIter], y[posIter], z[posIter]);
    if ((abs (radius-radMax) < 0.5) || (abs (radius-radMin) < 0.5))
      onSideSet[i] = true;
    
  }
  
  intensivePrint ("DONE GETTING SIDES");
  
}

void mesh::checkAndProject (std::vector<double> &v0, std::vector<double> &v1, 
                            std::vector<double> &v2, std::vector<double> &p0) {

  std::vector<double> orig (3, 0);
  std::vector<double> n0 = getNormalVector (v0, v1, v2);
  double dist            = projWonV_Dist (v0, n0, orig);
  double pRadius         = getRadius (p0[0], p0[1], p0[2]);  
  
  if (dist < 0) {
    n0[0] = n0[0] * (-1);
    n0[1] = n0[1] * (-1);
    n0[2] = n0[2] * (-1);
  }
  
  double tiny = 0.01;
  dist = projWonV_Dist (v0, n0, orig);
  // cout << "PROJ: "   << dist << endl;
  // cout << "BEFORE: " << getRadius (p0[0], p0[1], p0[2]) << endl;
  dist = projWonV_Dist (v0, n0, orig);
  if        (pRadius > dist && (abs (pRadius-radMax) < 1)) {
    double dif = abs (radMax - dist);
    if (abs(dif) > 1){
    cout << dist << " ABOVE " << pRadius << endl;
    cout << dif << " ABOVE" << endl;}
    p0[0] = p0[0] + (n0[0]) * (dif + tiny);
    p0[1] = p0[1] + (n0[1]) * (dif + tiny);
    p0[2] = p0[2] + (n0[2]) * (dif + tiny);    
  } else if (pRadius < dist && (abs (pRadius-radMin) < 1)) {
    double dif = abs (radMin - dist);
    cout << dist << " " << pRadius << endl;
    cout << dif << " BELOW" << endl;
    p0[0] = p0[0] + n0[0] * (dif+tiny) * (-1);
    p0[1] = p0[1] + n0[1] * (dif+tiny) * (-1);
    p0[2] = p0[2] + n0[2] * (dif+tiny) * (-1);        
  }
  // cout << "AFTER: " << getRadius (p0[0], p0[1], p0[2]) << endl;
  // cin.get();

  // ofstream myfile1 ("points.txt", ios::out);
  // myfile1 << p0[0] << ' ' << p0[1] << ' ' << p0[2] << endl;
  // myfile1.close ();
  //
  //
  // ofstream myfile2 ("facets.txt", ios_base::app);
  // myfile2 << v0[0] << ' ' << v0[1] << ' ' << v0[2] << endl;
  // myfile2 << v1[0] << ' ' << v1[1] << ' ' << v1[2] << endl;
  // myfile2 << v2[0] << ' ' << v2[1] << ' ' << v2[2] << endl;
                                       
}

double mesh::returnUpdate (vector<vector<double>> &vec, double &valMsh, 
                           size_t &reg, int &ind) {
                             
  // Small function that checks for the existance of a vector, and returns its value at [reg][pnt],
  // in addition to additing the value to vslMsh. If the vector does not exist, it just returns 
  // valMsh. This is the workhorse function tti model updates.
                             
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
                                   size_t &reg, int &ind) {
  
  // Small function that checks for the existance of a vector, and returns its value at [reg][pnt].
  // If it does not exist, it returns the default (in this case valmsh).
  
  if (not vec.empty ()) {
    return vec[reg][ind];
  } else {
    return valMsh;
  }
  
}

double mesh::returnUpdate1d (vector<vector<double>> &vec, double &valMsh, size_t &reg, int &ind, 
                             double &val1d) {
   
  // Overload that returns the parameter added to a 1d background model.
  
  if (not vec.empty ()) {
    return val1d + vec[reg][ind];
  } else {
    return valMsh;
  }
  
}               
                        
void mesh::dump (exodus_file &eFile) {
  
  eFile.writeVariable (c11, "du1");
  eFile.writeVariable (c66, "c66");
  
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
    
      // get 1d background model.
      double vs1d, vp1d, rho1d;
      background_models backgroundMod;
      if (mod.onedBackground == "europe_model")
        backgroundMod.eumod (rad, vs1d, vp1d, rho1d);
    
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
      rhoNew = returnUpdateAbsolute (mod.rho, rhoMsh, region, pnt);
      vsvNew = returnUpdateAbsolute (mod.vsv, vsvMsh, region, pnt);
      vshNew = returnUpdateAbsolute (mod.vsh, vshMsh, region, pnt);
      vpvNew = returnUpdateAbsolute (mod.vpv, vpvMsh, region, pnt);
      vphNew = returnUpdateAbsolute (mod.vph, vphMsh, region, pnt);
      
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