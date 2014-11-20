#include "classes.hpp"

using namespace std;

std::string model::returnDirection () {

/*
  Function returns the current direction of interpolation from the 
  parameter file.
*/

  return direction;

}

void model::findChunkCenters () {
  
/*
  Finds the center of model chunks.
*/
  
  vector<double> ctr (3, 0);
  xCtr.resize (worldSize);
  yCtr.resize (worldSize);
  zCtr.resize (worldSize);

  double xAvg=0;
  double yAvg=0;
  double zAvg=0;
  size_t totParam=0;
  for (size_t r=0; r<numModelRegions; r++) {

    size_t numParam = x[r].size ();
    for (size_t i=0; i<numParam; i++) {

      xAvg += x[r][i];
      yAvg += y[r][i];
      zAvg += z[r][i];
      totParam++;

    }
  }

  xCtr[myRank] = xAvg / totParam;
  yCtr[myRank] = yAvg / totParam;
  zCtr[myRank] = zAvg / totParam;
    
}

void model::findNeighbouringChunks () {
  
/*
  Finds neighbouring model chunks (for MPI). This assumes that the centers have already been
  calculated (using findChunkCenters ()). Use largeDistance to trick the algorithm against
  trying to call itself a neighbour.
*/
  
  size_t maxNeighbours = 8;
  size_t largeDistance = 1e10;
  
  sum1DVector (xCtr);
  sum1DVector (yCtr);
  sum1DVector (zCtr);
  
  vector<double> distArray (worldSize, 0);  
  for (size_t i=0; i<worldSize; i++) {
    
    double xDist = xCtr[myRank] - xCtr[i];
    double yDist = yCtr[myRank] - yCtr[i];
    double zDist = zCtr[myRank] - zCtr[i];
    
    if (getRadius (xDist, yDist, zDist) != 0) {
      distArray[i] = getRadius (xDist, yDist, zDist);                
    } else {
      distArray[i] = largeDistance;
    }
    
  }
  
  while (neighbourArray.size () < maxNeighbours) {
    
    size_t ind = getSmallestIndex (distArray);
    neighbourArray.push_back (ind);
    distArray[ind] = largeDistance;
    
    if (neighbourArray.size () == (worldSize - 1))
      break;
    
  }
  
  std::sort (neighbourArray.begin (), neighbourArray.end ());
  
}

void model::broadcastNeighbouringChunks () {
  
/*
  Function to broadcast all neighbouring chunk arrays, and append to initial model.
*/

  // Definitions.
  size_t numBroadcast  = vsh[0].size ();
  size_t numNeighbours = neighbourArray.size ();  
  int X_TAG            = 0;
  int Y_TAG            = 1;
  int Z_TAG            = 2;
  int KRN_TAG          = 3;
  originalSize         = numBroadcast;
  
  // Allocations.
  double *recX   = new double [numBroadcast];
  double *recY   = new double [numBroadcast];
  double *recZ   = new double [numBroadcast];
  double *recKrn = new double [numBroadcast];
  
  // For all neighbours previously found, broadcast full arrays to each processor.
  for (size_t i=0; i<numNeighbours; i++) {
  
    MPI::COMM_WORLD.Isend (&x[0][0],   numBroadcast, MPI::DOUBLE, neighbourArray[i], X_TAG);
    MPI::COMM_WORLD.Isend (&y[0][0],   numBroadcast, MPI::DOUBLE, neighbourArray[i], Y_TAG);
    MPI::COMM_WORLD.Isend (&z[0][0],   numBroadcast, MPI::DOUBLE, neighbourArray[i], Z_TAG);
    MPI::COMM_WORLD.Isend (&vsh[0][0], numBroadcast, MPI::DOUBLE, neighbourArray[i], KRN_TAG);
    
  }
  
  for (size_t i=0; i<numNeighbours; i++) {
    
    MPI::COMM_WORLD.Recv (recX,   numBroadcast, MPI::DOUBLE, neighbourArray[i], X_TAG);
    MPI::COMM_WORLD.Recv (recY,   numBroadcast, MPI::DOUBLE, neighbourArray[i], Y_TAG);
    MPI::COMM_WORLD.Recv (recZ,   numBroadcast, MPI::DOUBLE, neighbourArray[i], Z_TAG);
    MPI::COMM_WORLD.Recv (recKrn, numBroadcast, MPI::DOUBLE, neighbourArray[i], KRN_TAG);

    x[0].insert   (x[0].end (),   recX,   recX+numBroadcast);
    y[0].insert   (y[0].end (),   recY,   recY+numBroadcast);
    z[0].insert   (z[0].end (),   recZ,   recZ+numBroadcast);
    vsh[0].insert (vsh[0].end (), recKrn, recKrn+numBroadcast);
    
    MPI::COMM_WORLD.Barrier ();

  }
    
  delete [] recX;
  delete [] recY;
  delete [] recZ;
  delete [] recKrn;
  
}

void model::readParameterFile () {
  
  /*
    Read the default parameter file in ./mod/parameters.txt. TODO ignore comment lines.
  */
  
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
    
    if (paramName[i] == "overwrite_crust")
      overwriteCrust = paramValue[i];
    
    if (paramName[i] == "interpolate_all")
      interpolateAll = paramValue[i];
    
  }      
  
}

void model::resetParams () {
  
  vsh.clear ();
    
}

void model::allocateArrays () {

/*
  Allocates all the moduli arrays (for use with extraction).
*/
  
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
  
  /*
    Generate a kd-tree for the model volume.
  */
  
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
  
  /*
    If a rotation angle is specified, generate the rotation matrix, and apply it to all cartesian
  points in a model volume.
  */
  
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
    
    intensivePrint ("No rotation found.");
    
  }  
  
}


void model::findMinMaxRadius () {
  
  /*
    Find the minimum and maximum radius of a model volume. The function also determines which
  cem chunks to include. TODO move cemIncludeChunks to a seperate function.
  */
  
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

void model::findMinMaxCartesian () {
  
  /*
    Finds the minimum and maximum xyz points in a certain model volume.
  */
  
  xMin = x[0][0];
  xMax = x[0][0];
  yMin = y[0][0];
  yMax = y[0][0];
  zMin = z[0][0];
  zMax = z[0][0];
  
  for (size_t r=0; r<numModelRegions; r++) {
    
    size_t numParams = x[r].size ();
    for (size_t i=0; i<numParams; i++) {
      
      if (x[r][i] < xMin)
        xMin = x[r][i];
      
      if (x[r][i] > xMax)
        xMax = x[r][i];
      
      if (y[r][i] < yMin)
        yMin = y[r][i];
      
      if (y[r][i] > yMax)
        yMax = y[r][i];
      
      if (z[r][i] < zMin)
        zMin = z[r][i];
      
      if (z[r][i] > zMax)
        zMax = z[r][i];                        
      
    }
  }
  
}

bool model::checkBoundingBox (double &x, double &y, double &z) {
  
  /*
    Checks the (cartesian) bounds of a model volume, and return true if a requested point is inside,
  and false if it is outside.  
  */

  if (x <= xMax && x >= xMin &&
      y <= yMax && y >= yMin &&
      z <= zMax && z >= zMin) {
        
    return true;
    
  } else {
    
    return false;
    
  }
      
}

void model::construct () {
  
  /*
    Takes the elastic moduli and constructs a set of parameters from the moduli. 
  */
      
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
