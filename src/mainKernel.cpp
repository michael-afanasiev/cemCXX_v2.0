#include "classes.hpp"

using namespace std;

int main () {
  
  MPI::Init ();

  string exoFileName = "kernelMesh.ex2";

  specfem3d_globe modType;
  model *mod =& modType;
  
  exodus_file exo ("/Users/michaelafanasiev/Desktop/netcdfKernel/" + exoFileName,mod->regionNames);
  mesh msh (exo);
  msh.initializeKernel      (exo);
  msh.interpolateAndSmooth (*mod);

  if (MPI::COMM_WORLD.Get_rank () == 0) {
    exo.putVarParams ();
    exo.putVarNames ();
  }
  msh.dumpKernel (exo);

  msh.extract (*mod);
  mod->write ();

  // the mpi finalize seems to conflict with the mpi close in the hdf5 libraries. dumb.
  // replacing it here by an MPI_BARRIER instead.
  MPI::COMM_WORLD.Barrier ();

}