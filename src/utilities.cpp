#include "classes.hpp"

using namespace std;

rotation_matrix::rotation_matrix (float &a, float &x, float &y, float&z) {

  rot11 = cos(a) + (x * x) * (1 - cos(a));
  rot21 = z * sin(a) + x * y * (1 - cos(a));
  rot31 = y * sin(a) + x * z * (1 - cos(a));
  rot12 = x * y * (1 - cos(a)) - z * sin(a);
  rot22 = cos(a) + (y * y) * (1 - cos(a));
  rot32 = x * sin(a) + y * z * (1 - cos(a));
  rot13 = y * sin(a) + x * z * (1 - cos(a));
  rot23 = x * sin(a) + y * z * (1 - cos(a));
  rot33 = cos(a) + (z * x) * (1 - cos(a));

  rot23 = (-1) * rot23;
  rot31 = (-1) * rot31;

}

void rotation_matrix::rotate (float &xOld, float &yOld, float &zOld,
                              float &xNew, float &yNew, float &zNew) {
                           
  float xTmp = rot11 * xOld + rot21 * yOld + rot31 * zOld;
  float yTmp = rot12 * xOld + rot22 * yOld + rot32 * zOld;
  float zTmp = rot13 * xOld + rot23 * yOld + rot33 * zOld;
  
  xNew = xTmp;
  yNew = yTmp;
  zNew = zTmp;

}

float deg2Rad (float &deg) {
  
  float rad = deg * M_PI / 180.;
  return rad;
  
}

float rad2Deg (float &rad) {
  
  float deg = rad * 180. / M_PI;
  return deg;
  
}

// MPI Things.
int getRank () {
  
  return MPI::COMM_WORLD.Get_rank ();
  
}

void broadcastInteger (int &bInt) {
  
  MPI::COMM_WORLD.Bcast (&bInt, 1, MPI::INT, 0);
  
}

void broadcast1DVector (vector<int> &bVector) {
  
  bool broadcast=true;
  if (MPI::COMM_WORLD.Get_rank () == 0)
    if (bVector.empty ())
      broadcast = false;
      
  MPI::COMM_WORLD.Bcast (&broadcast, 1, MPI::BOOL, 0);
  
  if (broadcast) {
    
    int size=0;
    if (MPI::COMM_WORLD.Get_rank () == 0)
      size = bVector.size ();
    
    MPI::COMM_WORLD.Bcast (&size, 1, MPI::INT, 0);
    if (MPI::COMM_WORLD.Get_rank () >  0)
      bVector.resize (size);
    
    MPI::COMM_WORLD.Bcast (&bVector[0], size, MPI::INT, 0);
  
  }
  
}

void broadcast2DVector (vector<vector<float>> &bVector) {  
  vector<vector<float>>::iterator outer;
  bool broadcast=true;
  if (MPI::COMM_WORLD.Get_rank () == 0)
    if (bVector.empty ())
      broadcast = false;      
  MPI::COMM_WORLD.Bcast (&broadcast, 1, MPI::BOOL, 0);  
  if (broadcast) {
    vector<int> axesSize;
    int containerSize = bVector.size();
    MPI::COMM_WORLD.Bcast (&containerSize, 1, MPI::INT, 0);    
    if (MPI::COMM_WORLD.Get_rank () > 0)
      bVector.resize (containerSize);    
    if (MPI::COMM_WORLD.Get_rank () == 0) {
      for (outer=bVector.begin(); outer!=bVector.end(); ++outer) {          
        axesSize.push_back (outer->size());          
      }
    }    
    for (size_t i=0; i<containerSize; i++) {      
      int size=0;
      if (MPI::COMM_WORLD.Get_rank () == 0) {
        size = axesSize[i];
      }      
      MPI::COMM_WORLD.Bcast (&size, 1, MPI::INT, 0);      
      float *recvBuf = new float [size];
      if (MPI::COMM_WORLD.Get_rank () > 0)
        bVector[i].resize (size);
      if (MPI::COMM_WORLD.Get_rank () == 0)
        std::copy (bVector[i].begin (), bVector[i].end(), recvBuf);      
      MPI::COMM_WORLD.Bcast (&recvBuf[0], size, MPI::FLOAT, 0);      
      if (MPI::COMM_WORLD.Get_rank () > 0) 
        std::copy (recvBuf, recvBuf+size, bVector[i].begin ());      
      delete [] recvBuf;      
    }
  }  
}


// Messages.
void intensivePrint (string message) {
  
  if (MPI::COMM_WORLD.Get_rank () == 0)
    cout << yel << message << rst << endl;
  
}

void error (string message) {
  
  if (MPI::COMM_WORLD.Get_rank () == 0) 
    cout << red << "\n" << message << "\n" << rst << flush << endl;
  
  MPI::COMM_WORLD.Abort (EXIT_FAILURE);
  
}