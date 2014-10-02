#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "kdtree.h"
#include "mpi.h"

using namespace std;

// ###### misc. utilities. ######

// Colour codes for pretty writes.
static const char *rst = "\x1b[0m";
static const char *red = "\x1b[31m";
static const char *grn = "\x1b[32m";
static const char *yel = "\x1b[33m";
static const char *mgn = "\x1b[35m";
static const char *blu = "\x1b[36m";

// Message helper functions.
void error             (std::string);
void intensivePrint    (std::string);

// MPI helper functions.
void broadcast2DVector (std::vector<std::vector<float>>&);
void broadcast1DVector (std::vector<int>&);
void broadcastInteger  (int &);
int  getRank           ();

// Math functions
float deg2Rad (float &deg);
float rad2Deg (float &rad);
float getRadius (float &x, float &y, float &z);
float projWonV_Dist (float &x, float &y, float &z, std::vector<float> &v, std::vector<float> &x0);
std::vector<float> getNormalVector (std::vector<float> &A, std::vector<float> &B, 
                                    std::vector<float> &C);

// ###### global variables ######
const double R_EARTH = 6371.0;

class facet {
  
public:
  
    
  // nodes of polygon.
  vector<float> v0; 
  vector<float> v1;
  vector<float> v2;
  
  // normal.
  vector<float> n;
  
  // outside set.
  vector<int> region;
  vector<int> index;

  facet (vector<float> v0, vector<float> v1, vector<float> v2);
    
};

// ###### classes ######
class rotation_matrix;
class exodus_file;
class model;
class ses3d;


class model {

protected:
  
  int myRank;
  int worldSize;
  
  int numModelParams=0;
  int numModelRegions=0;
  
  std::vector<int> regionSize;
  
  // Spherical co-ordinate data.
  std::vector<std::vector<float>> col, lon, rad;
  
  // Cartesian co-ordinate data.
  std::vector<std::vector<float>> x, y, z;
    
  // Elastic moduli.
  std::vector<std::vector<float>> c11, c12, c13, c14, c15, c16;
  std::vector<std::vector<float>> c22, c23, c24, c25, c26, c33;
  std::vector<std::vector<float>> c34, c35, c36, c44, c45, c46;
  std::vector<std::vector<float>> c55, c56, c66;
  
  // density.
  vector<vector<float>> rho;

  // subset parameters.
  vector<vector<float>> vsv, vsh, vpv, vph;    
  vector<vector<float>> vsi, vpi;  
  
  // KD-trees.
  std::vector<std::vector<kdtree*>> trees;  
  std::vector<std::vector<int>> datKD;
  
  // Rotation parameters.
  float angle, xRot, yRot, zRot;
  
  // Model extremes.
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;
  float rMax, rMin;
  
  int xMinI[2], yMinI[2], zMinI[2];
  int xMaxI[2], yMaxI[2], zMaxI[2];
  
  
  std::vector<float> p1, p2, p3, p4;
  std::vector<float> n1, n2, n3, n4;
         
  virtual void read  (void) =0;
  virtual void write (void) =0;
  
  void createKDtree      ();
  void rotate            ();
  void findMinMax        ();
  void findEdgePlanes    ();
  void findConvexHull    ();
  void findMinMaxCartesian ();
  void findMinMaxRadius ();
  
  bool testEdge (float &x, float &y, float &z);

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
  
  void broadcast         ();
  void readFile          (std::vector<std::vector<float>> &vec, std::string type);
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
  int numNodes;
  int numElem;
  int numElemBlock;
  const int numNodePerElem=4;
  
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
  
  rotation_matrix (float &ang, float &x, float &y, float &z);
  void rotate     (float &x,    float &y,    float &z,
                   float &xNew, float &yNew, float &zNew);
  
private:
  
  float rot11, rot12, rot13, rot21, rot22, rot23;
  float rot31, rot32, rot33;

};