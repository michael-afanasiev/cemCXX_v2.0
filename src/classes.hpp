#include <iostream>

class exodus_file;

class exodus_file {
  
private:
  
  // define read/write characteristics for exodus files.
  int comp_ws = 8;
  int io_ws   = 0;
  int relert;
  int chrret;
  int ier;
  int idexo;
  
  // initialize with dummy filename for safety.
  std::string fileName;
  
public:
  
  exodus_file (std::string);

};