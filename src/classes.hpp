#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "kdtree.h"
#include "mpi.h"

using namespace std;

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
double deg2Rad                      (double &deg);
double rad2Deg                      (double &rad);
double getRadius                    (double &x, double &y, double &z);
double projWonV_Dist                (double &x, double &y, double &z, std::vector<double> &v, 
                                    std::vector<double> &x0);
                                    
void xyz2ColLonRad                 (double &x, double &y, double &z, double &col, double &lon, 
                                    double &rad);                                                                          
void colLonRad2xyz                 (double &x, double &y, double &z, double &col, double &lon, 
                                    double &rad);                                      
                                    
std::vector<double> getNormalVector (std::vector<double> &A, std::vector<double> &B, 
                                    std::vector<double> &C);

// ###### global variables ######
const double R_EARTH = 6371.0;

// ###### classes ######
class rotation_matrix;
class exodus_file;
class model;
class ses3d;


class model {

protected:
  
  int myRank;
  int worldSize;
  
  size_t numModelParams=0;
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
  std::vector<std::vector<kdtree*>> trees;  
  std::vector<std::vector<int>> datKD;
  
  // Rotation parameters.
  double angle, xRot, yRot, zRot;
  
  // Model extremes.
  double xMin, yMin, zMin;
  double xMax, yMax, zMax;
  double xCtr, yCtr, zCtr;
  double rMax, rMin;

  virtual void read  (void) =0;
  virtual void write (void) =0;
  
  void createKDtree        ();
  void rotate              ();
  void findMinMax          ();
  void findMinMaxCartesian ();
  void findMinMaxRadius    ();

};

class ses3d: public model {
  
public:
  
  ses3d (string path, string symSys);
  
protected:
    
  // book keeping.
  std::string path;
  std::string symSys;    
    
  void read  (void);
  void write (void) {};
  
  void readFile          (std::vector<std::vector<double>> &vec, std::string type);
  void broadcast         ();
  void convert2Cartesian ();
  void convert2Radians   ();
      
};

class exodus_file {
  
private:
  
  // define read/write characteristics for exodus files.
  float vers  = 0;
  int comp_ws = 8;
  int io_ws   = 0;
  int relert  = 0;
  int chrret  = 0;
  int idexo   = 0;
  int ier     = 0;
  
  // misc. mesh details.
  size_t numNodes;
  size_t numElem;
  size_t numElemBlock;
  const size_t numNodePerElem=4;
  
  // bookkeeping arrays.
  int *nodeNumMap;
  int *elemNumMap;
  int *blockNumMap;
  std::vector<int> connectivity;
  
  // initialize with dummy filename for safety.
  std::string fileName;
  
  // internal private functions.
  void exodusCheck      (int, std::string);
  void getInfo          ();
  void allocate         ();
  void getNodeNumMap    ();
  void getElemNumMap    ();
  void openFile         ();  
  void closeFile        ();
  void getConnectivity  ();
  
  int getNumElemInBlock (int &elmBlockId);
  
    
public:
  
  // Constructor.
  exodus_file   (std::string);
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