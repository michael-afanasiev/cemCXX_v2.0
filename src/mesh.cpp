#include "classes.hpp"

using namespace std;

mesh::mesh (exodus_file &eFile) {
  
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
    
    double rhoNew, vsvNew, vshNew, vpvNew, vphNew;
    
    // Initialize values to mesh values.
    double vshMsh = sqrt (c44[mshInd] / rho[mshInd]);
    double vsvMsh = sqrt (c55[mshInd] / rho[mshInd]);
    double vpvMsh = sqrt (c22[mshInd] / rho[mshInd]);
    double vphMsh = sqrt (c22[mshInd] / rho[mshInd]);
    double rhoMsh = rho[mshInd];  
  
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