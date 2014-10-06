#include "classes.hpp"

using namespace std;

edge::edge (vector<float> v0In, vector<float> v1In) {
  
  v0 = v0In;
  v1 = v1In;
  
}

facet::facet (vector<float> v0In, vector<float> v1In, vector<float> v2In, vector<float> ctr) {
  
  v0 = v0In;
  v1 = v1In;
  v2 = v2In;
  
  n = getNormalVector (v0, v1, v2);
  
  float zero=0.;
  // if ((projWonV_Dist (ctr[0], ctr[1], ctr[2], n, v0)) > 0) {
  if ((projWonV_Dist (zero, zero, zero, n, v0)) > 0) {
    n[0] = n[0] * -1;
    n[1] = n[1] * -1;
    n[2] = n[2] * -1;
  }
  
  outside = true;
  hasOutsideSet = false;
  remove = false;
    
}

void facet::reset () {
  
  oldRegionSet = regionSet;
  oldIndexSet  = indexSet;  
  dMax         = 0.; 
  
  regionSet.clear ();      
  indexSet.clear  ();
  
}

std::vector<edge> facet::findHorizon (facet &master) {

  std::vector<edge> freeEdge;
  
  // case: facet hangs on by two edge.
  if ((master.v0 == v0 || master.v0 == v1 || master.v0 == v2) && 
      (master.v1 == v0 || master.v1 == v1 || master.v1 == v2) &&
      (master.v2 == v0 || master.v2 == v1 || master.v2 == v2)) {
    
    float d0 = distFromLine (master.pMax, v0, v1);
    float d1 = distFromLine (master.pMax, v1, v2);
    float d2 = distFromLine (master.pMax, v2, v0);
    
    if (d0 > d1 && d0 > d2) 
      freeEdge.push_back (edge (v0, v1));
    
    if (d1 > d0 && d1 > d2)
      freeEdge.push_back (edge (v1, v2));
    
    if (d2 > d0 && d2 > d1)
      freeEdge.push_back (edge (v2, v0));
    
    return freeEdge;
  }
  
  // case: facet hangs on by one edge.
  if (master.v0 != v0 && master.v1 != v0 && master.v2 != v0) {
    freeEdge.push_back (edge (v0, v1));
    freeEdge.push_back (edge (v0, v2));
    return freeEdge;
  }

  if (master.v0 != v1 && master.v1 != v1 && master.v2 != v1) {
    freeEdge.push_back (edge (v0, v1));
    freeEdge.push_back (edge (v1, v2));
    return freeEdge;
  }

  if (master.v0 != v2 && master.v1 != v2 && master.v2 != v2) {
    freeEdge.push_back (edge (v0, v2));
    freeEdge.push_back (edge (v1, v2));
    return freeEdge;
  }
  
  // case: facet hangs on by two edges.
  
          
}

bool facet::checkNeighbour (facet &master) {
  
  vector<bool> verts (3,0);
  
  if (master.v0 == v0)
    verts[0] = true;
  
  if (master.v0 == v1)
    verts[0] = true;
  
  if (master.v0 == v2)
    verts[0] = true;
  
  if (master.v1 == v0)
    verts[1] = true;
  
  if (master.v1 == v1)
    verts[1] = true;  
    
  if (master.v1 == v2)
    verts[1] = true;
  
  if (master.v2 == v0)
    verts[2] = true;
  
  if (master.v2 == v1)
    verts[2] = true;
  
  if (master.v2 == v2)
    verts[2] = true;
  
  
  int nTrue = 0;
  for (size_t i=0; i<verts.size(); i++)
    if (verts[i])
      nTrue++;
  
  if (nTrue == 2 || nTrue == 3) {
    return true;
  } else {
    return false;
  }
  
  
}

pointTracker::pointTracker (model &mod) {
  
  pointTrack.resize (mod.numModelRegions);
  for (size_t r=0; r<mod.numModelRegions; r++) {
    
    int numParams = mod.x[r].size ();
    pointTrack[r].resize (numParams);
    for (size_t i=0; i<numParams; i++) {
      pointTrack[r][i] = true;                  
    }
  }
  
}

void pointTracker::setFalse (size_t &reg, size_t &ind) {
  
  pointTrack[reg][ind] = false;
  
}

bool pointTracker::check (size_t &reg, size_t &ind) {
  
  if (pointTrack[reg][ind] == true) {
    return true;
  } else {
    return false;
  }
  
}

void pointTracker::reset (model &mod) {
  
  for (size_t r=0; r<mod.numModelRegions; r++) {
    
    int numParams = mod.x[r].size ();
    for (size_t i=0; i<numParams; i++) {
      pointTrack[r][i] = true;                  
    }
  }
  
}