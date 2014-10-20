#include "classes.hpp"

int main () {
  
  MPI::Init ();

  ses3d modType;
  model *mod =& modType;
  
  discontinuity topo;
  
  std::vector<std::string> fileNames;
  std::vector<std::string>::iterator fileNameIter;
  
  fileNames = getRequiredChunks (*mod);
    
  for (fileNameIter=fileNames.begin (); fileNameIter!=fileNames.end(); ++fileNameIter) {
    
    exodus_file exo (*fileNameIter, mod->regionNames);  
  
    mesh msh (exo);
    
    if (mod->direction == "interpolate") {
      
      msh.interpolate (*mod);
      msh.dump (exo);      
      
    } else if (mod->direction == "extract") {
      
      msh.extract (*mod);
      
    } else if (mod->direction == "interpolate_topography") {
      
      msh.interpolateTopography (topo);
      msh.dump (exo);
      
    }
  }
  
  if (mod->direction == "extract")
    mod->write ();

  MPI::Finalize ();
    
}
