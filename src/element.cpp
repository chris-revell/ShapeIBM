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

element::element(const int& cell, const int& Label, const float& initialx, const float& initialy,const int& n1,const int& n2,const float& adhesion)  {
  pos(0)            =initialx;
  pos(1)            =initialy;
  initialpos(0)     =initialx;
  initialpos(1)     =initialy;
  label             = Label;
  parent            = cell;
  baselineadhesion = adhesion;
  //adhesionmagnitude = adhesion;
  neighbours.push_back(n1);
  neighbours.push_back(n2);
}

void element::SetAdhesion(void){
  adhesionmagnitude = baselineadhesion;
}

element::~element() {}
