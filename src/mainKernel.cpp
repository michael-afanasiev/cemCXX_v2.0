#include "classes.hpp"

using namespace std;

int main () {
  
  MPI::Init ();

  string exoFileName = "kernelMesh.ex2";

  specfem3d_globe modType;
  model *mod =& modType;
  
  // if (MPI::COMM_WORLD.Get_rank () == 0)
  if (MPI::COMM_WORLD.Get_rank () == 3) {
  exodus_file exo ("/Users/michaelafanasiev/Desktop/netcdfKernel/" + exoFileName, mod->regionNames);
  
  // MPI::COMM_WORLD.Barrier ();
  
  kernel kern (exo);
  kern.interpolate (*mod);

  exo.putVarParams ();
  exo.putVarNames ();
  
  kern.write (exo);
}
    
  MPI::Finalize ();

}