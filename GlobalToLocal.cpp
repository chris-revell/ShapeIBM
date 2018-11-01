//
//  GlobalToLocal.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 01/11/2018.
//
//

#include "cell.hpp"
#include "tissue.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

void GlobalToLocal(tissue& Tissue){
  for (int ii=0;ii<Tissue.Nc;ii++){
    Tissue.Cells[ii].xb = Tissue.xbglobal.cols(ii*Tissue.Nbcell,(ii+1)*Tissue.Nbcell-1);
  }
}
