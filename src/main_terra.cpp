#include "classes.hpp"

using namespace std;

int main (int argc, char *argv[]) {
  
  MPI::Init ();

  // specfem3d_globe modType;
  // ses3d modType;

  int rank_bump = atoi(argv[1]);
  int size=1024;
  if (MPI::COMM_WORLD.Get_rank() + rank_bump < size) {
      std::cout << MPI::COMM_WORLD.Get_rank() + rank_bump << std::endl;
  terragrid modType(rank_bump);
 
  model *mod =& modType;
  
  std::vector<std::string> fileNames;
  std::vector<std::string>::iterator fileNameIter;
  
  fileNames = getRequiredChunks (*mod);

  MPI::COMM_WORLD.Barrier ();

  if (MPI::COMM_WORLD.Get_rank () == 0)
    std::cout << grn << "Initialization complete." << rst << endl;

  for (fileNameIter=fileNames.begin (); fileNameIter!=fileNames.end(); ++fileNameIter) {
  
    exodus_file exo (*fileNameIter, mod->regionNames, mod->returnDirection ());    
    mesh msh            (exo, mod->returnDirection ());
    msh.initializeModel (exo);

    if (mod->direction == "interpolate") {
      
      msh.interpolate (*mod);
      msh.dump         (exo);      
      
    } else if (mod->direction == "extract") {
    
      msh.extract (*mod);
      
    } else if (mod->direction == "interpolate_topography") {

      discontinuity topo;
      msh.interpolateTopography (topo);
      msh.dump                   (exo);
      
    }

  }

  cout << "DONE! " << MPI::COMM_WORLD.Get_rank() << endl;
  MPI::COMM_WORLD.Barrier ();
  cout << "DONE! " << MPI::COMM_WORLD.Get_rank() << endl;
  if (mod->direction == "extract") mod->write ();
  }
 
  MPI::COMM_WORLD.Barrier ();
  MPI::Finalize ();
    
}
