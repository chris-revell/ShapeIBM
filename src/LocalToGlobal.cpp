//
//  LocalToGlobal.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#include "LocalToGlobal.hpp"
#include "element.hpp"
#include "cell.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

void LocalToGlobal(tissue& Tissue){
  int& Nb = Tissue.Nb;
  int& Nc = Tissue.Nc;
  vector<cell>& Cells = Tissue.Cells;

  Nb= 0;
  for (int ii=0;ii<Nc;ii++){
    Nb = Nb+Cells[ii].Nb;
  }
  if (Tissue.xbglobal.n_cols < Nb){
    Tissue.xbglobal.resize(2,Nb);
    Tissue.fbglobal.resize(2,Nb);
    Tissue.ubglobal.resize(2,Nb);
    Tissue.indices.resize(2,Nb);
  }
  for (int ii=0;ii<Nc;ii++){
    for (int jj=0; jj<Cells[ii].Nb;jj++){
      Tissue.xbglobal(0,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].pos(0);
      Tissue.xbglobal(1,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].pos(1);
      Tissue.fbglobal(0,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].fb(0);
      Tissue.fbglobal(1,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].fb(1);
      Tissue.indices(0,Cells[ii].Elements[jj].label)  = ii;
      Tissue.indices(1,Cells[ii].Elements[jj].label)  = jj;
    }
  }
  Tissue.ubglobal.zeros();
}
