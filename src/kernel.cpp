#include "classes.hpp"

using namespace std;

kernel::kernel (exodus_file &eFile) {
  
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

void kernel::interpolate (model &mod) {
  
  intensivePrint ("Interpolating.");
  size_t setSize = interpolatingSet.size ();
  int percent = (setSize / omp_get_max_threads ()) / 100.;
    
  int pCount = 0;
  int pIter  = 0;
  
  for (size_t i=0; i<setSize; i++) {

    // extract node number.
    size_t nodeNum = interpolatingSet[i] - 1;        
   
    // use du2 as a scratch array to avoid doubly visiting points.
    if (du1[nodeNum] == 1)
      continue;
         
    // find closest point [region specific].
    for (size_t r=0; r<mod.numModelRegions; r++) {
                                                
      kdres *set = kd_nearest3 (mod.trees[r], x[nodeNum], y[nodeNum], z[nodeNum]);
      void *ind  = kd_res_item_data (set);
      int point  = * (int *) ind;
      kd_res_free (set); 

      // get distance.
      double xDist = x[nodeNum] - mod.x[r][point];
      double yDist = y[nodeNum] - mod.y[r][point];
      double zDist = z[nodeNum] - mod.z[r][point];      
      double dist  = getRadius (xDist, yDist, zDist);

      // save distance.
      value[nodeNum] = dist;        
      
      // mark that we've visited here.
      du1[nodeNum] = 1;
      
    }
    
    if (omp_get_thread_num () == 0) {
      pCount++;
      if (pCount % percent == 0) {
        cout << pIter << " %\r" << flush;
        pIter++;
      }
    }
         
  }
  
  cout << grn << "Done." << rst << endl;
  
}

void kernel::write (exodus_file &eFile) {
  
  eFile.writeVariable (value, "krn");
  
  
}