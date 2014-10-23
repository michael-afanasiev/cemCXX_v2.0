#include "classes.hpp"
#include <netcdf>

using namespace std;

specfem3d_globe::specfem3d_globe () {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
  
  numModelRegions = 1;

  readParameterFile ();
  read              ();
  findMinMaxRadius  ();
  createKDtree      ();
  
}

void specfem3d_globe::read () {
  
  intensivePrint  ("Reading specfem3d_globe model.");
  readCoordNetcdf ();
  
  vph = readParamNetcdf ("alphahKernelCrustMantle.nc");  
  vpv = readParamNetcdf ("alphavKernelCrustMantle.nc");  
  vsh = readParamNetcdf ("betahKernelCrustMantle.nc");  
  vsv = readParamNetcdf ("betavKernelCrustMantle.nc");  
  eta = readParamNetcdf ("etaKernelCrustMantle.nc");
  rho = readParamNetcdf ("rhoKernelCrustMantle.nc");
  
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
    int numWroteProcs = procDim.getSize ();
    int numGLLPoints  = kernDim.getSize ();
  
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

void specfem3d_globe::readCoordNetcdf () {
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  
  // Filename.
  std::string fileName = path + "/xyzCrustMantle.nc";
  
  // Open file.
  try {
  
    NcFile dataFile (fileName, NcFile::read);
  
    // Get variable.
    NcVar NcRadius = dataFile.getVar ("radius");
    NcVar NcTheta  = dataFile.getVar ("theta");
    NcVar NcPhi    = dataFile.getVar ("phi");
  
    // Get array sizes.
    NcDim procDim     = NcRadius.getDim (0);
    NcDim coordDim    = NcRadius.getDim (1);
    int numWroteProcs = procDim.getSize ();
    size_t numGLLPoints  = coordDim.getSize ();
  
    if (myRank == 0)
      std::cout << mgn << "Number of solver processers:\t " << numWroteProcs
      << "\nNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;  
    
    if (numWroteProcs != worldSize) {
      error ("Wrong number of MPI processors.");
      exit (EXIT_FAILURE);
        
    }
  
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
      
    // Kernels are output as spherical coordinates.
    if (interpolationType == "kernel") {
    
      for (size_t i=0; i<numGLLPoints; i++) {

        double x, y, z;
        double col = yStore[i];
        double lon = zStore[i];
        double rad = xStore[i] * R_EARTH;
        colLonRad2xyz (x, y, z, col, lon, rad);
        xStore[i] = x;
        yStore[i] = y;
        zStore[i] = z;      
      
      }        
    
    }
    
    // copy to model arrays.
    x.resize (1); y.resize (1); z.resize (1);
    x[0].resize (numGLLPoints); y[0].resize (numGLLPoints); z[0].resize (numGLLPoints);
    std::copy (xStore, xStore+numGLLPoints, x[0].begin ());
    std::copy (yStore, yStore+numGLLPoints, y[0].begin ());
    std::copy (zStore, zStore+numGLLPoints, z[0].begin ());
    
    delete [] xStore;
    delete [] yStore;
    delete [] zStore;
    
  } catch (NcException &error) {

    std::cout << red << "Failure reading NetCDF file. See below" << endl;
    std::cout << error.what () << std::endl;         
    std::exit (EXIT_FAILURE);   

  }
  
}