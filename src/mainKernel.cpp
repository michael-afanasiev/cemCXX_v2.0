#include "classes.hpp"

using namespace std;

int main () {
  
  MPI::Init ();

//  string exoFileName = "kernelMesh.ex2";
 
  specfem3d_globe modType;
  model *mod =& modType;
  
  std::vector<std::string> fileNames;
  std::vector<std::string>::iterator fileNameIter;
  
  fileNames = getRequiredChunks (*mod);

  MPI::COMM_WORLD.Barrier ();
  
  for (fileNameIter=fileNames.begin (); fileNameIter!=fileNames.end(); fileNameIter++) {

    exodus_file exo (*fileNameIter, mod->regionNames, mod->returnDirection ());    
    mesh msh (exo, mod->returnDirection ());
    msh.initializeKernel      (exo);
    msh.findKernelInRange    (*mod);
    msh.interpolateAndSmooth (*mod);

    if (MPI::COMM_WORLD.Get_rank () == 0)
    exo.writeNew (*fileNameIter + ".kernel.ex2", msh, msh.krn);

  }

//  if (MPI::COMM_WORLD.Get_rank () == 0)
//    exo.writeNew ("./test.ex2", msh, msh.krn);

//  msh.extract (*mod);
// mod->write ();

  // the mpi finalize seems to conflict with the mpi close in the hdf5 libraries. dumb.
  // replacing it here by an MPI_BARRIER instead.
  MPI::COMM_WORLD.Barrier ();

}
