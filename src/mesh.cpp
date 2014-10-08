#include "classes.hpp"

using namespace std;

mesh::mesh (exodus_file &eFile) {
  
  eFile.getXYZ (x, y, z);
  c11 = eFile.getVariable ("du1");  
  connectivity = eFile.returnConnectivity ();
  
}

void mesh::interpolate (model &mod) {
  
  intensivePrint ("Interpolating.");
  size_t sizeConnect = connectivity.size ();
  std::vector<double> v0 (3, 0), v1 (3, 0), v2 (3, 0), v3 (3, 0);
  std::vector<double> c0 (3, 0), p0 (3, 0);
  int percent = sizeConnect / 100.;
  
  int pCount = 0;
  for (size_t i=0; i<sizeConnect; i+=numNodePerElem) {

    // get edges of tet.
    v0[0] = x[connectivity[i]];
    v0[1] = y[connectivity[i]];
    v0[2] = z[connectivity[i]];

    v1[0] = x[connectivity[i+1]];
    v1[1] = y[connectivity[i+1]];
    v1[2] = z[connectivity[i+1]];

    v2[0] = x[connectivity[i+2]];
    v2[1] = y[connectivity[i+2]];
    v2[2] = z[connectivity[i+2]];

    v3[0] = x[connectivity[i+3]];
    v3[1] = y[connectivity[i+3]];
    v3[2] = z[connectivity[i+3]];
    
    // get centroid of tet.
    c0[0] = (v0[0] + v1[0] + v2[0] + v3[0]) / 4.;
    c0[1] = (v0[1] + v1[1] + v2[1] + v3[1]) / 4.;
    c0[2] = (v0[2] + v1[2] + v2[2] + v3[2]) / 4.;

    // find closest point to centroid [region specific].
    for (size_t r=0; r<mod.numModelRegions; r++) {

      kdres *set = kd_nearest3 (mod.trees[r], c0[0], c0[1], c0[2]);
      void *ind  = kd_res_item_data (set);
      int point  = * (int *) ind;
      kd_res_free (set);

      p0[0] = mod.x[r][point];
      p0[1] = mod.y[r][point];
      p0[2] = mod.z[r][point];

      if (testInsideTet (v0, v1, v2, v3, p0)) {

        c11[connectivity[i]  -1] = 1;
        c11[connectivity[i+1]-1] = 1;
        c11[connectivity[i+2]-1] = 1;
        c11[connectivity[i+3]-1] = 1;        
        
        } else {

        c11[connectivity[i]  -1] = 0;
        c11[connectivity[i+1]-1] = 0;
        c11[connectivity[i+2]-1] = 0;
        c11[connectivity[i+3]-1] = 0;

      }
    }            
  }
  
}

void mesh::dump (exodus_file &eFile) {
  
  eFile.writeVariable (c11, "du1");
  
}