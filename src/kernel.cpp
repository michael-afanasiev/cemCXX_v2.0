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
  int percent = (setSize) / 100.;
  
  struct {
    double distance;
    int    rank;
  } distanceIn[numNodes], distanceOut[numNodes];
    
  int pCount       = 0;
  int pIter        = 0;  
  double searchRad = 100;  
  double *bufValue = new double [setSize];
  double *interpParam = new double [setSize];
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

      // cout << nodeNum << ' ' << setSize << endl;
      // Save distance and rank in struct.
      distanceIn[nodeNum].distance = dist;
      distanceIn[nodeNum].rank     = myRank;
      
      // save param.
      interpParam[nodeNum] = mod.vsh[r][point];

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
  
  // Figure out where the minimum distance is.
  MPI::COMM_WORLD.Allreduce (distanceIn, distanceOut, setSize, MPI_DOUBLE_INT, MPI_MINLOC);  
  
  for (size_t i=0; i<setSize; i++) {
    
    bufValue[i] = 0.;
    if (distanceOut[i].rank == myRank)      
      bufValue[i] = interpParam[i];              
    
  }
  
  if (myRank == 0) {
    MPI::COMM_WORLD.Reduce (MPI_IN_PLACE, bufValue, setSize, MPI_DOUBLE, MPI_SUM, 0);
  } else {
    MPI::COMM_WORLD.Reduce (bufValue,     bufValue, setSize, MPI_DOUBLE, MPI_SUM, 0);
  }
  
  std::copy (bufValue, bufValue+setSize, value.begin ());
  
  
}

void kernel::write (exodus_file &eFile) {
  
  if (myRank == 0)
    eFile.writeVariable (value, "krn");
  
  
}