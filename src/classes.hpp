#include <iostream>

class exodus_file;
class model_file;

// Global variables.
const double R_EARTH = 6371.0;

class model_file {

private:
  
  std::string modelFormat;

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
  
  // initialize with dummy filename for safety.
  std::string fileName;

  // internal private functions.
  void exodusCheck (int, std::string);
    
public:
  
  // Constructor.
  exodus_file (std::string);
  
  void openFile ();
  void closeFile ();
};
