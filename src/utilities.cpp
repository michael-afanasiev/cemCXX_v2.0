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

std::vector<float> getNormalVector (std::vector<float> &A,
                                    std::vector<float> &B,
                                    std::vector<float> &C) {
                                      
   // Gets the normal vector to 3 points in 3-dimensions. Used to determine the equation of a plane.
                                      
  std::vector<float> AB;
  std::vector<float> AC;
  std::vector<float> n;
  
  AB.resize (3);
  AC.resize (3);
  n.resize  (3);
  
  AB[0] = B[0] - A[0];
  AB[1] = B[1] - A[1];
  AB[2] = B[2] - A[2];
  
  AC[0] = C[0] - A[0];
  AC[1] = C[1] - A[1];
  AC[2] = C[2] - A[2];
  
  n[0] = AB[1] * AC[2] - AB[2] * AC[1];
  n[1] = AB[0] * AC[2] - AB[2] * AC[0] * (-1);
  n[2] = AB[0] * AC[1] - AB[1] * AC[0];
  
  float magnitude = sqrt (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] = n[0] / magnitude;
  n[1] = n[1] / magnitude;
  n[2] = n[2] / magnitude;
  
  return n;
  
}

float projWonV_Dist (float &x, float &y, float&z, std::vector<float> &v, std::vector<float> &x0) {
  
  // Projects a vector x - x0 onto the plane v.
  
  float dotVW = v[0] * (x-x0[0]) + v[1] * (y-x0[1]) + v[2] * (z - x0[2]);
  float magV  = sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  
  return dotVW / magV;
  
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

float getRadius (float &x, float &y, float &z) {
  
  float rad = sqrt (x*x + y*y + z*z);
  return rad;
  
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