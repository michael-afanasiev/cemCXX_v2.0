#include "classes.hpp"

using namespace std;

kernel::kernel (exodus_file &eFile, model &mod) {
  
  myRank    = MPI::COMM_WORLD.Get_rank ();
  worldSize = MPI::COMM_WORLD.Get_size ();
  
  eFile.getXYZ (x, y, z);
  
  connectivity     = eFile.returnConnectivity ();
  nodeNumMap       = eFile.returnNodeNumMap   ();
  numNodes         = eFile.numNodes;
  interpolatingSet = eFile.returnInterpolatingSet ();
  
  du1.resize (numNodes);
  value.resize (numNodes);
  
}

void kernel::write (exodus_file &eFile) {
  
  if (myRank == 0)
    eFile.writeVariable (value, "krn");
    
}