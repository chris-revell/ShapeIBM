//
//  element.cpp
//  ImmersedBoundary
//
//  Created by <author> on 07/11/2018.
//
//

#include "element.hpp"

element::element(const float& initialx, const float& initialy)  {
  pos(0)=initialx;
  pos(1)=initialy;
  neighbours[0]=n1;
  neighbours[1]=n2;
}


element::~element() {}
