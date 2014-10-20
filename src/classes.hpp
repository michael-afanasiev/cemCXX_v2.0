#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>

#include <sys/stat.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <omp.h>

#include "kdtree.h"
#include "mpi.h"

using namespace std;

// ###### classes ######
class background_models;
class rotation_matrix;
class discontinuity;
class elasticTensor;
class exodus_file;
class model;
class ses3d;
class mesh;


// ###### misc. utilities. ######

// Colour codes for pretty writes.
const std::string rst = "\x1b[0m";
const std::string red = "\x1b[31m";
const std::string grn = "\x1b[32m";
const std::string yel = "\x1b[33m";
const std::string mgn = "\x1b[35m";
const std::string blu = "\x1b[36m";

// Message helper functions.
void error             (std::string);
void intensivePrint    (std::string);

// MPI helper functions.
void broadcast2DVector (std::vector<std::vector<double>>&);
void broadcast1DVector (std::vector<int>&);
void broadcastInteger  (size_t &);
void broadcastInteger  (int &);
int  getRank           ();

// Math functions
void xyz2ColLonRad     (double &x, double &y, double &z, double &col, double &lon, double &rad);                                                                         
void colLonRad2xyz     (double &x, double &y, double &z, double &col, double &lon, double &rad);

bool testInsideTet     (vector<double> &v0, vector<double> &v1, vector<double> &v2, 
  vector<double> &v3, vector<double> &p0,double &l1, double &l2, double &l3, double &l4);
                       
size_t vectorSize2d    (std::vector<std::vector<double>>);                        

double deg2Rad         (double &deg);
double rad2Deg         (double &rad);
double getRadius       (double &x, double &y, double &z);
double projWonV_Dist   (std::vector<double> &x, std::vector<double> &v, std::vector<double> &x0);
double interpolateTet  (std::vector<double> &vec, size_t &n0, size_t &n1, size_t &n2, size_t &n3, 
  double &l0, double &l1, double &l2, double &l3); 
                       
std::vector<double> getNormalVector   (std::vector<double> &A, std::vector<double> &B, std::vector<double> &C);
std::vector<double> returnVector      (double &x, double &y, double &z);     
std::vector<string> getRequiredChunks (model &mod);    

// ###### global variables ######
const double R_EARTH = 6371.0;
const double CLOSE   = 1;
const double TINY    = 0.01;
const double BIGTINY = 0.1;
const double RAD_400 = 5971.0;
const double RAD_670 = 5701.0;
const double RAD_CMB = 3480.0;
const double RAD_ICB = 1221.0;

class model {
  
  friend class mesh;

protected:
  
  int myRank;
  int worldSize;
  
  // Parameter file.
  std::string path;
  std::string symSys;      
  std::string interpolationType;
  std::string convert_to_1_second;
  std::string onedBackground;
 
  // Number of model regions.
  size_t numModelRegions=0;
  std::vector<int> regionSize;
  
  // Spherical co-ordinate data.
  std::vector<std::vector<double>> col, lon, rad;
  
  // Cartesian co-ordinate data.
  std::vector<std::vector<double>> x, y, z;
    
  // Elastic moduli.
  std::vector<std::vector<double>> c11, c12, c13, c14, c15, c16;
  std::vector<std::vector<double>> c22, c23, c24, c25, c26, c33;
  std::vector<std::vector<double>> c34, c35, c36, c44, c45, c46;
  std::vector<std::vector<double>> c55, c56, c66;
  
  // density.
  vector<vector<double>> rho;

  // subset parameters.
  vector<vector<double>> vsv, vsh, vpv, vph;    
  vector<vector<double>> vsi, vpi;  
  
  // KD-trees.
  std::vector<kdtree*> trees;  
  std::vector<std::vector<int>> datKD;
  
  // Rotation parameters.
  double angle, xRot, yRot, zRot;
  
  // Model extremes [physical].
  std::vector<double> xMin, yMin, zMin;
  std::vector<double> xMax, yMax, zMax;
  std::vector<double> minRadRegion, maxRadRegion;

  // Center of box.
  double xCtr, yCtr, zCtr;  
 

  // Internal functions.
  void rotate            ();
  void construct         ();
  void findMinMax        ();
  void createKDtree      ();
  void allocateArrays    ();
  void findBoundingBox   ();
  void findMinMaxRadius  ();  
  void readParameterFile ();
  
  int testBoundingBox (double x, double y, double z);
  
public:
  

  virtual void read  (void) =0;
  virtual void write (void) =0;
  
  std::vector<std::string> regionNames;
  std::set<std::string> colChunks;
  std::set<std::string> lonChunks;
  std::vector<double> rMax, rMin;  
  std::string meshDirectory;
  std::string direction;    

};

class attenuation {

protected:
  
  int    nRelaxationMechanisms=3;
  double tau_s[3];
  double D[3];

public:
  
  double QL6 (double &rad);
  double correctQL6 (double &rad);

};

class ses3d: public model {
  
public:
  
  ses3d ();
    
  void read  (void);
  void write (void);
  
protected:
  
  void readFile          (std::vector<std::vector<double>> &vec, std::string type);
  void writeFile         (std::vector<std::vector<double>> &vec, std::string type);
  void broadcast         ();
  void convert2Cartesian ();
  void convert2Radians   ();
      
};


class mesh {

  friend class exodus_file;
  friend class model;
  
public:
  
  mesh (exodus_file &);  
  
  void dump        (exodus_file &);
  void extract     (model &);
  void interpolate (model &);
  void interpolateTopography (discontinuity &topo);

protected:
  
  // co-ordinates.
  std::vector<double> x, y, z;  
 
  // exodus file name
  std::string eFileName;
  
  // extremes
  double xMin, xMax, yMin, yMax, zMin, zMax;
  double radMin, radMax;
  
  // Elastic moduli.
  std::vector<double> c11, c12, c13, c14, c15, c16;
  std::vector<double> c22, c23, c24, c25, c26, c33;
  std::vector<double> c34, c35, c36, c44, c45, c46;
  std::vector<double> c55, c56, c66, elv;
  
  // Density
  std::vector<double> rho;
  
  // connectivity array.
  std::vector<int> connectivity;
  
  // node number map.
  std::vector<int> nodeNumMap;
  
  // bool array for region finding.
  std::vector<int> interpolatingSet;
  
  // vector for side sets.
  std::vector<int> sideSetSide;
  std::vector<int> sideSetElem;
    
  // kdtree.
  kdtree *tree;
  std::vector<int> datKD;
    
  // Keep track of side sets.
  std::vector<bool> onSideSet;
    
  // misc. mesh details.
  const size_t numNodePerElem=4;
  size_t numNodes;
  
  elasticTensor breakdown     (model &mod, double &x, double &y, double &z, 
                               size_t &region, size_t &mshInd, int &point);                                
  double returnUpdateAbsolute (vector<vector<double>> &vec, double &valMsh, 
                               size_t &reg, int &pnt);
  double returnUpdate1d       (vector<vector<double>> &vec, double &valMsh, 
                               size_t &reg, int &pnt, double &val1d);
  double returnUpdate         (vector<vector<double>> &vec, double &valMsh, 
                               size_t &reg, int &pnt);
  double SBTRKTUpdate         (vector<vector<double>> &vec, double &valMsh, 
                               size_t &reg, int &pnt);
  void createKDTree ();
  void getMinMaxDimensions ();
  bool checkBoundingBox (double &x, double &y, double &z);
  bool checkInterpolatingRegion (double &x, double &y, double &z, double minRad, double maxRad);
  void getSideSets ();
  
  void checkAndProject (std::vector<double> &v0, std::vector<double> &v1,
                        std::vector<double> &v2, std::vector<double> &p0);
  
};

class elasticTensor {

public:
  
  double c11, c12, c13, c14, c15, c16;
  double c22, c23, c24, c25, c26, c33;
  double c34, c35, c36, c44, c45, c46;
  double c55, c56, c66, rho;

};

class exodus_file {
  
  friend class mesh;
  
protected:
  
  // define read/write characteristics for exodus files.
  float vers  = 0;
  int comp_ws = 8;
  int io_ws   = 0;
  int relert  = 0;
  int chrret  = 0;
  int idexo   = 0;
  int ier     = 0;
  
  // misc. mesh details.
  size_t numElem;
  size_t numNodes;
  size_t numNodeSets;
  size_t numSideSets;
  size_t numElemBlock;
  const size_t numNodePerElem=4;
  const size_t numSidePerElem=3;
  
  // bookkeeping arrays.
  std::vector<int> nodeNumMap;
  int *elemNumMap;
  int *blockNumMap;
  int *nodeSetNumMap;
  int *sideSetNumMap;
  std::vector<int>  connectivity;
  std::vector<int>  interpolatingSet;
  std::vector<int>  sideSetSide;
  std::vector<int>  sideSetElem;
  std::vector<bool> onSideSet;
  
  // initialize with dummy filename for safety.
  std::string fileName;
    
  // internal private functions.
  void getInfo          ();
  void allocate         ();
  void getNodeNumMap    ();
  void getElemNumMap    ();
  void openFile         ();  
  void closeFile        ();
  void getSideSets      ();
  void getConnectivity  (std::vector<std::string> regionNames);
  void getNodeSets      (std::vector<std::string> regionNames);
  void writeVariable    (std::vector<double> &var, std::string varName);
  void exodusCheck      (int, std::string);
  void getXYZ           (std::vector<double> &x, std::vector<double> &y, 
    std::vector<double> &z);
      
  int getNumElemInBlock (int &elmBlockId);
  int getNumNodeInSet   (int &nodeSetId);

  std::vector<bool>   returnOnSideSet        ();
  std::vector<int>    returnConnectivity     ();
  std::vector<int>    returnNodeNumMap       ();
  std::vector<int>    returnInterpolatingSet ();
  std::vector<int>    returnSideSetSide      ();
  std::vector<int>    returnSideSetElem      ();
  std::vector<double> getVariable            (std::string varName);

  std::string returnName ();
  

      
public:
  
  // Constructor.
  exodus_file   (std::string, std::vector<std::string>);
  ~exodus_file  ();
  
  void printMeshInfo  ();
  
};


class rotation_matrix {

public:
  
  rotation_matrix (double &ang, double &x, double &y, double &z);
  void rotate     (double &x,    double &y,    double &z,
                   double &xNew, double &yNew, double &zNew);
  
private:
  
  double rot11, rot12, rot13, rot21, rot22, rot23;
  double rot31, rot32, rot33;

};

class background_models {

public:
  
  void eumod                    (double &, double &, double &, double &);
  void prem_no220               (double &, double &, double &, double &);
  void eumod_vpPrem_vsPremLt670 (double &, double &, double &, double &);
          
};

class discontinuity {

public:
  
  discontinuity ();
  
  std::vector<double> elv;  
  std::vector<int> KDdat;  
  kdtree *tree;
  
  void readTopography ();

};
