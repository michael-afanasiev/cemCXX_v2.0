#include <iostream>
#include <vector>
#include <fstream>
#include "mpi.h"

// Colour codes for pretty writes.
static const char *rst = "\x1b[0m";
static const char *red = "\x1b[31m";
static const char *grn = "\x1b[32m";
static const char *yel = "\x1b[33m";
static const char *mgn = "\x1b[35m";
static const char *blu = "\x1b[36m";
// End colour codes.

class exodus_file;
class model_file;

// Global variables.
const double R_EARTH = 6371.0;

class model_file {
  
public:
  
  model_file (std::string, std::string, std::string);

private:
  
  std::string modFormat;
  std::string symSys;
  std::string modPath;
  
  std::vector<std::vector<float>> colattitude;
  std::vector<std::vector<float>> longitude;
  std::vector<std::vector<float>> radius;
  
  // TTI model parameters.
  std::vector<std::vector<float>> vsh;
  std::vector<std::vector<float>> vsv;
  std::vector<std::vector<float>> vph;
  std::vector<std::vector<float>> vpv;
  std::vector<std::vector<float>> rho;
       
  void readSES3D (std::vector<std::vector<float>>&, std::string);
  void broadcastModelDefinitions ();
  void broadcastTTI ();

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
  const int numNodePerElem = 4;
  
  // bookkeeping arrays.
  int *nodeNumMap;
  int *elemNumMap;
  int *blockNumMap;
  int *connectivity; //or
  std::vector<std::vector<int>> connectivityVec;
  
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
  
    
public:
  
  // Constructor.
  exodus_file   (std::string);
  ~exodus_file  ();
  
  void printMeshInfo  ();
  
};
