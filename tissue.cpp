//
//  tissue.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "cell.hpp"

using namespace arma;
using namespace std;

tissue::tissue(const int& GridSize,const int& dimensions) {
  Ng=GridSize;
  xg = cube(Ng+1,Ng+1,2,fill::zeros);// Fluid grid array
  sg = mat(Ng+1,Ng+1,fill::zeros);
  xbglobal = mat(2,0,fill::zeros);   // Positions of all boundary points in the system
  std::vector<cell> Cells;           // Vector of cell objects in system
  Nc = 0;                            // Number of cells in the system
  xmin=-static_cast<float>(dimensions);                  // Fluid domain dimensions
  xmax=static_cast<float>(dimensions);                   // Fluid domain dimensions
  hg=float(xmax-xmin)/float(Ng);     // Fluid mesh width

  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
      for (int jj=0;jj<Ng+1;jj++){
          xg(ii,jj,0)=xmin+ii*hg;
          xg(ii,jj,1)=xmin+jj*hg;
      }
  }
}

void tissue::AddCell(const int& BoundPoints, const float& len, const float& initialx, const float& initialy) {
  Nb = BoundPoints;
  Cells.push_back(cell(Nb,len,initialx,initialy));
  Nc++;
}

/*
void tissue::CombineBoundaryPositions(){
  if (xbglobal.n_cols()<Nc*Nb){
    xbglobal.resize(2,Nc*Nb
  }
  for (int ii=0;ii<Nc;ii++){
    xbglobal.rows(ii*Nb,(ii+1)*Nb-1) = Cells[ii].xb
  }
}*/

//sources;  // Positions of all fluid sources in the system
//sinks;    // Positions of all fluid sinks in the system


tissue::~tissue() {}
