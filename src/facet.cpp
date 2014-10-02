#include "classes.hpp"

using namespace std;

facet::facet (vector<float> v0In, vector<float> v1In, vector<float> v2In) {
  
  v0 = v0In;
  v1 = v1In;
  v2 = v2In;
  
  n = getNormalVector (v0, v1, v2);
  
}