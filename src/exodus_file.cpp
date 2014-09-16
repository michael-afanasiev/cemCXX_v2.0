#include <exodusII.h>
#include "classes.hpp"

// Class constructor for exodus file.
exodus_file::exodus_file (std::string fname) {
  
  fileName = fname;
  
}

// Opens an exodus file and populates the idexo field.
void exodus_file::openFile () {
  
  std::cout << "Opening exodus file: " << fileName << std::flush << std::endl;
  idexo = ex_open (fileName.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
  
}

// Closes the exodus file.
void exodus_file::closeFile () {
  
  exodusCheck (ex_close (idexo), "ex_close");
  std::cout << "File closed succesfully." << std::flush << std::endl;
  
}

// Checks to make sure we're doing a sane operation on the exodus file.
void exodus_file::exodusCheck (int ier, std::string function) {
  
  if (ier != 0) {
    std::cout << "Exodus library error in " << function << std::flush 
      << std::endl;
    exit (EXIT_FAILURE);
  }
  
}