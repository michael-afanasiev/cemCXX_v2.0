#include "classes.hpp"

int main () {
  
  MPI::Init ();

  ses3d modType ("./mod", "tti_noRho_isoVp");
  model *mod =& modType;
  
  if (MPI::COMM_WORLD.Get_rank () == 0)
  exodus_file exo ("col000-090.lon000-090.rad6361-6371.000.ex2");
  
  
  MPI::Finalize ();
  
}