#include "classes.hpp"

using namespace std;

void model::readParameterFile () {
  
  ifstream inputFile ("./mod/parameters.txt");
  vector<string> paramName, paramValue;
  string paramNameDum, paramValueDum;
  
  while (inputFile >> paramNameDum >> paramValueDum) {
    
    paramName.push_back  (paramNameDum);
    paramValue.push_back (paramValueDum);
    
  }
  
  for (size_t i=0; i<paramName.size (); i++) {
    
    if (paramName[i] == "path")
      path = paramValue[i];
    
    if (paramName[i] == "symmetry_system")
      symSys = paramValue[i];
    
    if (paramName[i] == "interpolation_type")
      interpolationType = paramValue[i];
    
    if (paramName[i] == "convert_to_1_second")
      convert_to_1_second = paramValue[i];
    
    if (paramName[i] == "one_d_background")
      onedBackground = paramValue[i];
    
    if (paramName[i] == "mesh_directory")
      meshDirectory = paramValue[i];
    
    if (paramName[i] == "direction")
      direction = paramValue[i];
    
    if (paramName[i] == "interpolating_region")
      regionNames.push_back (paramValue[i]);
    
    if (paramName[i] == "taper")
      taper = paramValue[i];
    
  }      
  
}

void model::resetParams () {
  
  vsh.clear ();
    
}

void model::allocateArrays () {

  // Allocates all the moduli arrays (for use with extraction).
  
  if (interpolationType != "kernel") {
    
    c11.resize (numModelRegions);
    c12.resize (numModelRegions);
    c13.resize (numModelRegions);
    c14.resize (numModelRegions);
    c15.resize (numModelRegions);
    c16.resize (numModelRegions);
    c22.resize (numModelRegions);
    c23.resize (numModelRegions);
    c24.resize (numModelRegions);
    c25.resize (numModelRegions);
    c26.resize (numModelRegions);
    c33.resize (numModelRegions);
    c34.resize (numModelRegions);
    c35.resize (numModelRegions);
    c36.resize (numModelRegions);
    c44.resize (numModelRegions);
    c45.resize (numModelRegions);
    c46.resize (numModelRegions);
    c55.resize (numModelRegions);
    c56.resize (numModelRegions);
    c66.resize (numModelRegions);
    rho.resize (numModelRegions);
    
  } else if (interpolationType == "kernel") {
    
    krn.resize (numModelRegions);
    
  }
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size ();
    
    if (interpolationType != "kernel") {
      
      c11[r].resize (numParams);
      c12[r].resize (numParams);
      c13[r].resize (numParams);
      c14[r].resize (numParams);
      c15[r].resize (numParams);
      c16[r].resize (numParams);
      c22[r].resize (numParams);
      c23[r].resize (numParams);
      c24[r].resize (numParams);
      c25[r].resize (numParams);
      c26[r].resize (numParams);
      c33[r].resize (numParams);
      c34[r].resize (numParams);
      c35[r].resize (numParams);
      c36[r].resize (numParams);
      c44[r].resize (numParams);
      c45[r].resize (numParams);
      c46[r].resize (numParams);
      c55[r].resize (numParams);
      c56[r].resize (numParams);
      c66[r].resize (numParams);
      rho[r].resize (numParams);  
      
    } else if (interpolationType == "kernel") {
      
      krn[r].resize (numParams);
      
    }
    
  }
    
}

void model::createKDtree () {
  
  intensivePrint ("Creating KD-trees.");
  
  trees.reserve (numModelRegions);
  datKD.resize  (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    kdtree *tree  = kd_create (3);
    size_t numParams = x[r].size ();
    datKD[r].resize (numParams);

#pragma omp parallel for firstprivate (r)
    for (size_t i=0; i<numParams; i++) {
      
      datKD[r][i] = i;
      kd_insert3 (tree, x[r][i], y[r][i], z[r][i], &datKD[r][i]);    
      
    }
    
    trees.push_back (tree);
  
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
  
  double ZERO           = 0.0;
  double NINETY         = 90.0;
  double NEG_NINETY     = -90.0;
  double ONE_EIGHTY     = 180.0;
  double NEG_ONE_EIGHTY = -180.0;
  
  rMin.resize (numModelRegions);
  rMax.resize (numModelRegions);
  
  for (size_t r=0; r<numModelRegions; r++) {
    rMin[r] = getRadius (x[r][0], y[r][0], z[r][0]);
    rMax[r] = getRadius (x[r][0], y[r][0], z[r][0]);
  }
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size();
    for (size_t i=0; i<numParams; i++) {
      
      double colDum, lonDum, radDum;
      xyz2ColLonRad (x[r][i], y[r][i], z[r][i], colDum, lonDum, radDum);
      
      if (radDum < rMin[r])
        rMin[r] = radDum;
      
      if (radDum > rMax[r])
        rMax[r] = radDum;
      
      if (colDum <= deg2Rad (NINETY))
        colChunks.insert ("col000-090.");
      
      if (colDum >= deg2Rad (NINETY))
        colChunks.insert ("col090-180.");
    
      if (lonDum >= deg2Rad (ZERO) && lonDum <= deg2Rad (NINETY))
        lonChunks.insert ("lon000-090.");

      if (lonDum >= deg2Rad (NINETY) && lonDum <= deg2Rad (ONE_EIGHTY))
        lonChunks.insert ("lon090-180.");

      if (lonDum <= deg2Rad (ZERO) && lonDum >= deg2Rad (NEG_NINETY))
        lonChunks.insert ("lon270-360.");

      if (lonDum <= deg2Rad (NEG_NINETY) && lonDum >= deg2Rad (NEG_ONE_EIGHTY))
        lonChunks.insert ("lon180-270.");
      
    }
  }
    
}

void model::construct () {
      
  if (symSys.compare (0, 3, "tti") == 0) {
        
    vsh.resize (numModelRegions);
    vsv.resize (numModelRegions);
    vph.resize (numModelRegions);
    vpv.resize (numModelRegions);

    for (size_t r=0; r<numModelRegions; r++) {

      size_t numParams = x[r].size ();
      vsh[r].resize (numParams);
      vsv[r].resize (numParams);
      vph[r].resize (numParams);
      vpv[r].resize (numParams);

    }    
  }

  for (size_t r=0; r<numModelRegions; r++) {

    size_t numParams = x[r].size ();
    for (size_t i=0; i<numParams; i++) {

      if (symSys.compare (0, 3, "tti") == 0) {                

        vsh[r][i] = sqrt (c44[r][i] / rho[r][i]);
        vsv[r][i] = sqrt (c55[r][i] / rho[r][i]);
        vpv[r][i] = sqrt (c22[r][i] / rho[r][i]);
        vph[r][i] = sqrt (c22[r][i] / rho[r][i]);

      }                  
    }
  }  
}
