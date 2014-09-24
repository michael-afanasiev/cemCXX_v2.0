#include "classes.hpp"

int main () {
  
  MPI::Init ();
  
  std::cout << rst << "Hello world." << std::endl;
  model_file SES3DMODEL ("SES3D", "TTI", 
    "/Users/michaelafanasiev/Development/src/code/comprehensive_earth_model/Europe");
  exodus_file exoFile ("col000-090.lon000-090.rad6361-6371.000.ex2");
  exoFile.printMeshInfo ();
  
  MPI::Finalize ();
  
}