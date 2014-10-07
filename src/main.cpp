#include "classes.hpp"

int main () {
  
  MPI::Init ();

  ses3d modType ("./mod", "tti_noRho_isoVp");
  model *mod =& modType;
  
  exodus_file exo ("col000-090.lon090-180.rad6361-6371.000.ex2");  
  
  mesh msh (exo);
  msh.interpolate (*mod);
  msh.dump (exo);
  
  MPI::Finalize ();
    
}