//
//  element.cpp
//  ImmersedBoundary
//
//  Created by <author> on 07/11/2018.
//
//

#include "element.hpp"
#include <iostream>
using namespace std;
using namespace arma;

element::element(const float& v0, const float& v1,const int& cell, const int& Totalb, const float& initialx, const float& initialy,const int& n1,const int& n2)  {
  pos           = vec(2,fill::zeros);
  internalforce = vec(2,fill::zeros);
  ub            = vec(2,fill::zeros);
  pos(0)        =initialx;
  pos(1)        =initialy;
  ub(0)         = v0;
  ub(1)         = v1;
  //internalforce(0)=0;
  //internalforce(1)=0;
  neighbours.push_back(n1);
  neighbours.push_back(n2);
  label         = Totalb;
  parent        = cell;
}


element::~element() {}
