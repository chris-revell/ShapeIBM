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
  int celllabel,elementlabel;
  for (int ii=0;ii<Tissue.Nb;ii++){
    celllabel    = Tissue.indices(0,ii);
    elementlabel = Tissue.indices(1,ii);
    Tissue.Cells[celllabel].Elements[elementlabel].pos = Tissue.xbglobal.col(ii);
    //Tissue.Cells[cellindex].Elements[elementindex].pos(1) = Tissue.xbglobal(1,ii);
    //Tissue.Cells[cellindex].Elements[elementindex].ub(0)  = Tissue.ubglobal(0,ii);
    //Tissue.Cells[cellindex].Elements[elementindex].ub(1)  = Tissue.ubglobal(1,ii);
  }
}
