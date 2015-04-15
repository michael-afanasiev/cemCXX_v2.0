#include "classes.hpp"
#include <netcdf>
#include <netcdf.h>
#include <netcdf_par.h>

using namespace std;

specfem3d_globe::specfem3d_globe () {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
  

  readParameterFile ();
  
  if (interpolationType == "kernel") {
    numModelRegions = 1;
  } else {
    numModelRegions = 3;
  }
  
  read                        ();
  adjustRegions               ();
  findMinMaxRadius            ();
  findMinMaxCartesian         ();
  findChunkCenters            ();

  if (interpolationType == "kernel") {
    findNeighbouringChunks      ();
    broadcastNeighbouringChunks ();
    createKDtree                ();
  }

  allocateArrays              ();
  
}

void specfem3d_globe::adjustRegions () {

  for (size_t r=0; r<numModelRegions; r++) {

    size_t numGLLpoints = x[r].size ();
    for (size_t i=0; i<numGLLpoints; i++) {
    
      double colLoc, lonLoc, radLoc;
      xyz2ColLonRad (x[r][i], y[r][i], z[r][i], colLoc, lonLoc, radLoc);
      
      if (r == 0) {

        if (radLoc <= RAD_CMB) 
          radLoc = RAD_CMB + TINY;
        
        if (radLoc > R_EARTH)
          radLoc = R_EARTH - TINY;

        if (abs (radLoc - RAD_400) < CLOSE)
          radLoc = RAD_400 - CLOSE;

        if (abs (radLoc - RAD_670) < CLOSE)
          radLoc = RAD_670 - CLOSE;

      } else if (r == 1) {

        if (radLoc >= RAD_CMB)
          radLoc = RAD_CMB - TINY;

        if (radLoc <= RAD_ICB)
          radLoc = RAD_ICB + TINY;

      } else if (r == 2) {

        if (radLoc > 1115) 
          radLoc = 1115 - TINY;

        if (radLoc < 5) {
          radLoc = 5 + TINY;
        }
      }

//      if (colLoc == 0.)
//        colLoc = deg2Rad(*TINY);
//
//      if (colLoc == 180.)
//        colLoc = 180 - deg2Rad(TINY);
//
//      if (lonLoc == 0.)
//        lonLoc = deg2Rad(TINY);
        
      colLonRad2xyz (x[r][i], y[r][i], z[r][i], colLoc, lonLoc, radLoc);

    }
    
  }

}

void specfem3d_globe::read () {
  
  intensivePrint  ("Reading specfem3d_globe model.");
  readCoordNetcdf ("xyz_reg01.nc");
  MPI::COMM_WORLD.Barrier ();
  
  if (interpolationType != "kernel") {
    readCoordNetcdf ("xyz_reg02.nc");
    MPI::COMM_WORLD.Barrier ();
    readCoordNetcdf ("xyz_reg03.nc");
    MPI::COMM_WORLD.Barrier ();
  }

  if (interpolationType == "kernel") {

//    vsh = readParamNetcdf ("alphahKernelCrustMantle.nc");  
    vsh = readParamNetcdf ("alphavKernelCrustMantle.nc");  
//    vsh = readParamNetcdf ("betahKernelCrustMantle.nc");
//    vsv = readParamNetcdf ("betavKernelCrustMantle.nc");  
//    eta = readParamNetcdf ("etaKernelCrustMantle.nc");
//    rho = readParamNetcdf ("rhoKernelCrustMantle.nc");

  }
  
}

void specfem3d_globe::write () {
  
  if (interpolationType == "kernel") {

    writeParamNetcdf (krn[0], "test.nc");

  } else {

    construct ();

    stringstream myRankStringStream;
    myRankStringStream << std::setw(6) << std::setfill ('0') << myRank;
    string myRankString = myRankStringStream.str ();

    for (size_t r=0; r<numModelRegions; r++) {
      for (size_t i=0; i<rho[r].size (); i++) {

        if (rho[r][i] == 0) {

          cout << "X Y Z: " << x[r][i] << ' ' << y[r][i] << ' ' << z[r][i] << endl;
          cout << "RHO IS ZERO AT: " << getRadius (x[r][i], y[r][i], z[r][i]);

        }
      }
    }

    writeParamNetcdfSerial (vph[0], "vph_reg01.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vpv[0], "vpv_reg01.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsh[0], "vsh_reg01.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsv[0], "vsv_reg01.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (rho[0], "rho_reg01.proc" + myRankString + ".nc");

    MPI::COMM_WORLD.Barrier ();
    
    writeParamNetcdfSerial (vph[1], "vph_reg02.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vpv[1], "vpv_reg02.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsh[1], "vsh_reg02.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsv[1], "vsv_reg02.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (rho[1], "rho_reg02.proc" + myRankString + ".nc");
    
    MPI::COMM_WORLD.Barrier ();
    writeParamNetcdfSerial (vph[2], "vph_reg03.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vpv[2], "vpv_reg03.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsh[2], "vsh_reg03.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (vsv[2], "vsv_reg03.proc" + myRankString + ".nc");
    writeParamNetcdfSerial (rho[2], "rho_reg03.proc" + myRankString + ".nc");

    MPI::COMM_WORLD.Barrier ();
  }

}

void specfem3d_globe::writeParamNetcdfSerial (vector<double> &vec, std::string fileName) {

  using namespace netCDF;
  using namespace netCDF::exceptions;

  // compression filters (9 is highest)
  bool enableShuffleFilter = true;
  bool enableDeflateFilter = true;
  int deflateLevel = 9;

  fileName = path + fileName;
  int totSize = vec.size ();
  double *writeArr = &vec[0];

  try {

    NcFile output (fileName, NcFile::replace);
    NcDim dDim = output.addDim ("param", totSize);

    vector <NcDim> dims;
    dims.push_back (dDim);
    NcVar data = output.addVar ("data", ncDouble, dims);

    data.setCompression (enableShuffleFilter, enableDeflateFilter, deflateLevel);

    data.putVar (writeArr);

    delete [] writeArr;

  } catch (NcException &error) {

    std::cout << error.what() << std::endl;
    std::cout << red << "Failure writing: " << fileName << std::endl;
    std::exit (EXIT_FAILURE);    

  }
}

void specfem3d_globe::writeParamNetcdf (vector<double> &vec, std::string fileName) {
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  
  int ncid, kernDimId, procDimId, varId, ids[2];
  size_t start[2], count[2];
  
  size_t numParam  = vec.size ();  
  double *writeArr = &vec[0];

  fileName = path + fileName;
  fileSavePrint (fileName);

  try {
   
    string checkstr = "./test.nc";
    nc_create     (fileName.c_str (), NC_NETCDF4|NC_MPIIO, &ncid);// MPI::COMM_WORLD, MPI::INFO_NULL, &ncid); 
    nc_def_dim    (ncid, "glob", numParam,  &kernDimId);
    nc_def_dim    (ncid, "proc", worldSize, &procDimId);
    
    ids[0] = procDimId;
    ids[1] = kernDimId; 
   
    if (interpolationType == "kernel") {
      nc_def_var (ncid, "rawKernel", NC_DOUBLE, 2, ids, &varId);
   } else {
      nc_def_var (ncid, "param", NC_DOUBLE, 2, ids, &varId);
    }

    nc_enddef  (ncid);
    
    start[0] = myRank;
    start[1] = 0;
    count[0] = 1;
    count[1] = numParam;
    
    nc_put_vara_double (ncid, varId, start, count, writeArr);
    
    nc_close (ncid);
  } catch (NcException &error) {    
    
    std::cout << error.what() << std::endl;
    std::cout << red << "Failure writing: " << fileName << std::endl;
    std::exit (EXIT_FAILURE);    
    
  }
  
}

vector<vector<double>> specfem3d_globe::readParamNetcdf (std::string fName) {
  
  // Open the kernel file output from the solver.
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  
  string fileName = path + fName;
  
  if (myRank == 0)
    std::cout << "Opening kernel file: " << blu << fileName << rst << " with " << worldSize 
      << " processors." << std::flush << std::endl;
  
  double *rawKernel;
  vector<vector<double>> vec;
  
  try {

    // Open the file.
    NcFile dataFile (fileName, NcFile::read);
    
    // Get variable.
    NcVar NcKernel = dataFile.getVar ("rawKernel");
    
    // Get array sizes.
    NcDim procDim = NcKernel.getDim (0);
    NcDim kernDim = NcKernel.getDim (1);  
    size_t numWroteProcs = procDim.getSize ();
    size_t numGLLPoints  = kernDim.getSize ();
  
    if (myRank == 0)
      std::cout << mgn << "Number of solver processers:\t " << numWroteProcs
        << "\nNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;
    
    // Set up the MPI read chunk array.
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.resize(2);
    count.resize(2);
    
    // Row major MPI read. Start at [myRank, 0]
    start[0] = myRank;
    start[1] = 0;
    
    // Read until end of line [myrank, numGLLPoints]
    count[0] = 1;
    count[1] = numGLLPoints;
    
    // Of course only read in with the number of processors used to create the file.
    if (myRank < numWroteProcs) {
      rawKernel = new double [numGLLPoints];
      NcKernel.getVar (start, count, rawKernel);
    }
      
    vec.resize (1);
    vec[0].resize (numGLLPoints);
  
    std::copy (rawKernel, rawKernel+numGLLPoints, vec[0].begin ());    
    return vec;      
        
    // Destructor will close file.        
    
  } catch (NcException &error) {
    
    std::cout << error.what() << std::endl;
    std::cout << red << "Failure reading: " << fileName << std::endl;
    std::exit (EXIT_FAILURE);
    
  }
    
}

void specfem3d_globe::readCoordNetcdf (std::string fileName) {
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  

  intensivePrint ("Reading: " + fileName);

  // Filename.
  fileName = path + fileName;
  
  // Open file.
  try {
  
    NcFile dataFile (fileName, NcFile::read);
  
    // Get variable.

    NcVar NcRadius; NcVar NcTheta; NcVar NcPhi;
    if (interpolationType == "kernel_old") {
      
      NcRadius = dataFile.getVar ("radius");      
      NcTheta  = dataFile.getVar ("theta");
      NcPhi    = dataFile.getVar ("phi");

    } else {

      NcRadius = dataFile.getVar ("x");
      NcTheta  = dataFile.getVar ("y");
      NcPhi    = dataFile.getVar ("z");
      
    }
  
    // Get array sizes.
    NcDim procDim     = NcRadius.getDim (0);
    NcDim coordDim    = NcRadius.getDim (1);
    size_t numWroteProcs = procDim.getSize ();
    size_t numGLLPoints  = coordDim.getSize ();
  
    if (myRank == 0)
      std::cout << mgn << "Number of solver processers:\t " << numWroteProcs
      << "\nNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;  
    
//    if (numWroteProcs != worldSize) {
//      error ("Wrong number of MPI processors.");
//      exit (EXIT_FAILURE);
        
//    }
  
    // Set up the MPI read chunk array.
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.resize(2);
    count.resize(2);
  
    // Row major MPI read. Start at [myRank, 0]
    start[0] = myRank;
    start[1] = 0;
  
    // Read until end of line [myrank, numGLLPoints]
    count[0] = 1;
    count[1] = numGLLPoints;
  
    // Preallocate cartesian arrays.
    double *xStore = new double [numGLLPoints];
    double *yStore = new double [numGLLPoints];
    double *zStore = new double [numGLLPoints];
  
    // Of course only read in with the number of processors used to create the file.
    
    if (myRank < numWroteProcs) {
        
      NcRadius.getVar (start, count, xStore);
      NcTheta.getVar  (start, count, yStore);
      NcPhi.getVar    (start, count, zStore);
    
    } 
      
    // Kernels are output as spherical coordinates, while actual model parameters are output as
    // normalized cartesian co-ordinates (why?). This takes care of the conversion.
    for (size_t i=0; i<numGLLPoints; i++) {
      if (interpolationType == "kernel_old") {

        double x, y, z;
        double col = yStore[i];
        double lon = zStore[i];
        double rad = xStore[i] * R_EARTH;
        colLonRad2xyz (x, y, z, col, lon, rad);
        xStore[i] = x;
        yStore[i] = y;
        zStore[i] = z;      
      
      } else {
        
        xStore[i] = xStore[i];
        yStore[i] = yStore[i];
        zStore[i] = zStore[i];
        
      }    
    }
    
    // copy to model arrays.
    std::vector<double> xDum, yDum, zDum;
    xDum.resize (numGLLPoints); yDum.resize (numGLLPoints); zDum.resize (numGLLPoints);
    std::copy (xStore, xStore+numGLLPoints, xDum.begin ());
    std::copy (yStore, yStore+numGLLPoints, yDum.begin ());
    std::copy (zStore, zStore+numGLLPoints, zDum.begin ());
    
    delete [] xStore;
    delete [] yStore;
    delete [] zStore;

    x.push_back (xDum);
    y.push_back (yDum);
    z.push_back (zDum);
    
  } catch (NcException &error) {

    std::cout << red << "Failure reading NetCDF file. See below" << endl;
    std::cout << error.what () << std::endl;         
    std::exit (EXIT_FAILURE);   

  }
  
  if (fileName.find ("reg01") != string::npos) {
    maxRadRegion.push_back (R_EARTH);
    minRadRegion.push_back (RAD_CMB);

#ifdef VERBOSE
    cout << "DONE reg01 " << myRank << endl;
#endif
  } else if (fileName.find ("reg02") != string::npos) {
    maxRadRegion.push_back (RAD_CMB);
    minRadRegion.push_back (RAD_ICB);
#ifdef VERBOSE
    cout << "DONE reg02 " << myRank << endl;
#endif
  } else if (fileName.find ("reg03") != string::npos) {
    maxRadRegion.push_back (RAD_ICB);
    minRadRegion.push_back (0.);
#ifdef VERBOSE
    cout << "DONE reg03 " << myRank << endl;
#endif
  }

}
