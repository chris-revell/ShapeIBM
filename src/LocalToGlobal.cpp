//
//  LocalToGlobal.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#include "LocalToGlobal.hpp"
#include "element.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

void LocalToGlobal(tissue& Tissue){

  for (int ii=0;ii<Tissue.Nb;ii++){

      Tissue.xbglobal.col(ii) = Tissue.Elements[ii].pos;
      Tissue.fbglobal.col(ii) = Tissue.Elements[ii].fb;      
    }
  Tissue.ubglobal.zeros();        
}
