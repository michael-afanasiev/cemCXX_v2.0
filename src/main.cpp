#include "classes.hpp"

int main () {
  
  std::cout << "Hello world." << std::endl;
  exodus_file exoFile ("col000-090.lon000-090.rad6361-6371.000.ex2");
  exoFile.printMeshInfo ();
  
}