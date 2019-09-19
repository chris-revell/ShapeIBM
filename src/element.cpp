//
//  element.cpp
//  ImmersedBoundary
//
//  Created by <author> on 07/11/2018.
//
//

#include "element.hpp"
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

element::element(const int& Label, const float& initialx, const float& initialy, const float& accum, const float& r0, const float& r1){
  pos(0)              = initialx;
  pos(1)              = initialy;
  fb.zeros();
  accumulatedEffector = accum;
  initialpos(0)       = initialx;
  initialpos(1)       = initialy;
  re0                 = r0;
  re1                 = r1;
}

element::~element() {}
