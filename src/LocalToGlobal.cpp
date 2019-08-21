//
//  LocalToGlobal.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#include "LocalToGlobal.hpp"
#include "element.hpp"
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

void LocalToGlobal(const vector<element>& Elements,mat& xbglobal,mat& fbglobal,mat& ubglobal,const int& Nb){

  for (int ii=0;ii<Nb;ii++){
      xbglobal.col(ii) = Elements[ii].pos;
      fbglobal.col(ii) = Elements[ii].fb;
    }
  ubglobal.zeros();
}
