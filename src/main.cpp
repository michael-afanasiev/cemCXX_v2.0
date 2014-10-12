#include "classes.hpp"

int main () {
  
  MPI::Init ();

  ses3d modType;
  model *mod =& modType;
  
  std::vector<std::string> fileNames;
  std::vector<std::string>::iterator fileNameIter;
  
  fileNames = getRequiredChunks (*mod);
    
  for (fileNameIter=fileNames.begin (); fileNameIter!=fileNames.end(); ++fileNameIter) {
    
    cout << *fileNameIter << endl;
    exodus_file exo (*fileNameIter);  
  
    mesh msh (exo);
    // msh.interpolate (*mod);
    msh.extract (*mod);
    msh.dump (exo);
    
  }
  
  MPI::Finalize ();
    
}
