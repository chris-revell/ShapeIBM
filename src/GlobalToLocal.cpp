//
//  GlobalToLocal.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 01/11/2018.
//
//

#include "tissue.hpp"
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

void GlobalToLocal(vector<element>& Elements,const mat& xbglobal,const int& Nb){
  for (int ii=0;ii<Nb;ii++){
    Elements[ii].pos = xbglobal.col(ii);
  }
}
