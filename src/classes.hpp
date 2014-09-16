#include <iostream>
#include <vector>
#include <fstream>

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
  int *connectivity;
  
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
