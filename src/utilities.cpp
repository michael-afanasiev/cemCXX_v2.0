#include "classes.hpp"

using namespace std;

rotation_matrix::rotation_matrix (double &a, double &x, double &y, double&z) {

  rot11 = cos(a) + (x * x) * (1 - cos(a));
  rot21 = z * sin(a) + x * y * (1 - cos(a));
  rot31 = (-1) * y * sin(a) + x * z * (1 - cos(a));
  rot12 = x * y * (1 - cos(a)) - z * sin(a);
  rot22 = cos(a) + (y * y) * (1 - cos(a));
  rot32 = x * sin(a) + y * z * (1 - cos(a));
  rot13 = y * sin(a) + x * z * (1 - cos(a));
  rot23 = (-1) * x * sin(a) + y * z * (1 - cos(a));
  rot33 = cos(a) + (z * z) * (1 - cos(a));

}

void rotation_matrix::rotate (double &xOld, double &yOld, double &zOld,
                              double &xNew, double &yNew, double &zNew) {
                           
  double xTmp = rot11 * xOld + rot21 * yOld + rot31 * zOld;
  double yTmp = rot12 * xOld + rot22 * yOld + rot32 * zOld;
  double zTmp = rot13 * xOld + rot23 * yOld + rot33 * zOld;
  
  xNew = xTmp;
  yNew = yTmp;
  zNew = zTmp;

}

std::vector<string> getRequiredChunks (model &mod) {
  
  std::set<std::string>::iterator colIter;
  std::set<std::string>::iterator lonIter;
  
  std::vector<std::string> modelChunks;
  
  if (mod.direction == "interpolate_topography") {
    
    mod.colChunks.clear (); mod.lonChunks.clear (); mod.rMin.clear ();
    
    mod.colChunks.insert ("col000-090.");
    mod.colChunks.insert ("col090-180.");
    
    mod.lonChunks.insert ("lon000-090.");
    mod.lonChunks.insert ("lon090-180.");
    mod.lonChunks.insert ("lon180-270.");
    mod.lonChunks.insert ("lon270-360.");
    
    mod.rMin.push_back (R_EARTH-100.);
    
  }
  
  DIR *dp = opendir (mod.meshDirectory.c_str ());
  struct dirent *dirp;
  if (dp == NULL)
    error ("Directory where you told me to find the exodus files doesn't exist.");
    
  while ((dirp = readdir (dp))) {
        
    std::string fileName = dirp->d_name;
    if ( fileName.find_last_of ("ex2") != std::string::npos) {

      for (colIter=mod.colChunks.begin (); colIter!=mod.colChunks.end (); ++colIter) {
        for (lonIter=mod.lonChunks.begin (); lonIter!=mod.lonChunks.end (); ++lonIter) {

          double radMax = stod (fileName.substr (30, 4));
                        
          if ((radMax >= *std::min_element (mod.rMin.begin (), mod.rMin.end ())) && 
              (fileName.find (*colIter) != std::string::npos) &&
              (fileName.find (*lonIter) != std::string::npos)) {

            modelChunks.push_back (mod.meshDirectory + "/" + fileName);

          }
                            
        }
      }
    }                
  }
  
  return modelChunks;
  
}

bool testInsideTet (vector<double> &v0, 
                    vector<double> &v1, 
                    vector<double> &v2, 
                    vector<double> &v3,
                    vector<double> &p0,
                    double &l1, double &l2, double &l3, double &l4) {
                      
                      
  double x1 = v0[0]; double y1 = v0[1]; double z1 = v0[2];
  double x2 = v1[0]; double y2 = v1[1]; double z2 = v1[2];
  double x3 = v2[0]; double y3 = v2[1]; double z3 = v2[2];
  double x4 = v3[0]; double y4 = v3[1]; double z4 = v3[2];
  double xp = p0[0]; double yp = p0[1]; double zp = p0[2];
                      
  double vecX = xp - x4;
  double vecY = yp - y4;
  double vecZ = zp - z4;

  double a = x1 - x4;
  double d = y1 - y4;
  double g = z1 - z4;
  double b = x2 - x4;
  double e = y2 - y4;
  double h = z2 - z4;
  double c = x3 - x4;
  double f = y3 - y4;
  double i = z3 - z4;
  
  double det = 1 / (( a * ( e * i - f * h ) ) - ( b * ( i * d - f * g ) ) +
    ( c * ( d * h - e * g ) ));
  
  double ai = det * (e * i - f * h);
  double bi = det * (d * i - f * g) * (-1);
  double ci = det * (d * h - e * g);
  double di = det * (b * i - c * h) * (-1);
  double ei = det * (a * i - c * g);
  double fi = det * (a * h - b * g) * (-1);
  double gi = det * (b * f - c * e);
  double hi = det * (a * f - c * d) * (-1);
  double ii = det * (a * e - b * d);
  
  l1 = ai * vecX + di * vecY + gi * vecZ;
  l2 = bi * vecX + ei * vecY + hi * vecZ;
  l3 = ci * vecX + fi * vecY + ii * vecZ;
  l4 = 1 - l1 - l2 - l3;
  
  if (l1 >= 0 && l2 >= 0 && l3 >= 0 && l4 >= 0) {
    return true;
  } else {
    return false;
  }
                      
}

std::vector<double> getNormalVector (std::vector<double> &A,
                                     std::vector<double> &B,
                                     std::vector<double> &C) {
                                      
   // Gets the normal vector to 3 points in 3-dimensions.
                                      
  std::vector<double> AB;
  std::vector<double> AC;
  std::vector<double> n;
  
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
  n[1] = (AB[0] * AC[2] - AB[2] * AC[0]) * (-1);
  n[2] = AB[0] * AC[1] - AB[1] * AC[0];
  
  double magnitude = sqrt (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] = n[0] / magnitude;
  n[1] = n[1] / magnitude;
  n[2] = n[2] / magnitude;
  
  return n;
  
}

double interpolateTet (std::vector<double> &vec, size_t &n0,  size_t &n1, size_t &n2, size_t &n3,
                       double &l0, double &l1, double &l2, double &l3) {
                         
  return l0 * vec[n0] + l1 * vec[n1] + l2 * vec[n2] + l3 * vec[n3];                         
                         
}

void xyz2ColLonRad (double &x, double &y, double &z, double &col, double &lon, double &rad) {

  rad = getRadius (x, y, z);
  col = acos (z / rad);
  lon = atan2 (y, x);

}



void colLonRad2xyz (double &x, double &y, double &z, double &col, double &lon, double &rad) {
  
  x = rad * cos (lon) * sin (col);
  y = rad * sin (lon) * sin (col);
  z = rad * cos (col);
  
}

double projWonV_Dist (std::vector<double> &x, std::vector<double> &v, std::vector<double> &x0) {
  
  // Projects a vector x - x0 onto the plane v.
  
  double dotVW = v[0] * (x0[0]-x[0]) + v[1] * (x0[1]-x[1]) + v[2] * (x0[2] - x[2]);
  
  return dotVW;
  
}

double deg2Rad (double &deg) {
  
  double rad = deg * M_PI / 180.;
  return rad;
  
}

double rad2Deg (double &rad) {
  
  double deg = rad * 180. / M_PI;
  return deg;
  
}

// MPI Things.
int getRank () {
  
  return MPI::COMM_WORLD.Get_rank ();
  
}

void broadcastString (std::string message) {

  cout <<"HI" << endl;
  int strLength;
  if (MPI::COMM_WORLD.Get_rank () == 0) 
    strLength = message.length ();
 
  MPI::COMM_WORLD.Bcast (&strLength, 1, MPI::INT, 0);

  const char *buf = message.c_str ();
  MPI::COMM_WORLD.Bcast (&buf, strLength, MPI::CHAR, 0);
  
}


void broadcastInteger (size_t &bInt) {
  
  // initialize variable.
  if (MPI::COMM_WORLD.Get_rank () > 0) 
    bInt = 1;
  
  MPI::COMM_WORLD.Bcast (&bInt, 1, MPI::UNSIGNED, 0);
  
}

void broadcastInteger (int &bInt) {
  
  // initialize variable.
  if (MPI::COMM_WORLD.Get_rank () > 0) 
    bInt = 1;
  
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

void broadcast1DVector (vector<double> &bVector) {
  
  bool broadcast=true;
  if (MPI::COMM_WORLD.Get_rank () == 0)
    if (bVector.empty ())
      broadcast = false;
      
  MPI::COMM_WORLD.Bcast (&broadcast, 1, MPI::BOOL, 0);
  
  if (broadcast) {
    
    int size=0;
    if (MPI::COMM_WORLD.Get_rank () == 0)
      size = bVector.size ();
    
    MPI::COMM_WORLD.Bcast (&size, 1, MPI::DOUBLE, 0);
    if (MPI::COMM_WORLD.Get_rank () >  0)
      bVector.resize (size);
    
    MPI::COMM_WORLD.Bcast (&bVector[0], size, MPI::DOUBLE, 0);
  
  }
  
}

void broadcast2DVector (vector<vector<double>> &bVector) {  
  vector<vector<double>>::iterator outer;
  bool broadcast=true;
  if (MPI::COMM_WORLD.Get_rank () == 0)
    if (bVector.empty ())
      broadcast = false; 
             
  MPI::COMM_WORLD.Bcast (&broadcast, 1, MPI::BOOL, 0);  
  if (broadcast) {
    vector<int> axesSize;
    size_t containerSize = bVector.size();
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
      double *recvBuf = new double [size];
      if (MPI::COMM_WORLD.Get_rank () > 0)
        bVector[i].resize (size);
      if (MPI::COMM_WORLD.Get_rank () == 0)
        std::copy (bVector[i].begin (), bVector[i].end(), recvBuf);      
      MPI::COMM_WORLD.Bcast (&recvBuf[0], size, MPI::DOUBLE, 0);      
      if (MPI::COMM_WORLD.Get_rank () > 0) 
        std::copy (recvBuf, recvBuf+size, bVector[i].begin ());      
      delete [] recvBuf;      
    }
  }
    
}

size_t vectorSize2d (std::vector<std::vector<double>> vec) {
  
  size_t totalSize =0;
  for (size_t r=0; r<vec.size (); r++) {
    
    totalSize += vec[r].size ();

  }
  
  return totalSize;
  
}

std::vector<double> returnVector (double &x, double &y, double &z) {
  
  // Helper function to return an allocated vector to a point in 3d space.
  std::vector<double> vec (3, 0);
  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
  
  return vec;    
  
}

double getRadius (double &x, double &y, double &z) {
  
  double rad = sqrt (x*x + y*y + z*z);
  
  // if (rad > 6371)
  //   rad = 6371;
  
  return rad;
  
}

// Messages.
void intensivePrint (string message) {
  
  if (MPI::COMM_WORLD.Get_rank () == 0)
    cout << yel << message << rst << endl;
  
}

void fileSavePrint (string message) {

  if (MPI::COMM_WORLD.Get_rank () == 0)
    cout << "Saving: " << blu << message << rst << endl;

}


void error (string message) {
  
  if (MPI::COMM_WORLD.Get_rank () == 0) 
    cout << red << "\n" << message << "\n" << rst << flush << endl;
  
  MPI::COMM_WORLD.Abort (EXIT_FAILURE);
  
}
