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

using namespace arma;
using namespace std;

tissue::tissue(const int& GridSize,const int& dimensions,const float& sourcestrength,const float& density,const float& viscosity,const float& timestep,const float& cen,const float& adhesion,const float& tension,const float& cellheight,const float& re){
  Ng       = GridSize;
  //Nbcell   = boundarypoints;
  Nb       = 0;
  dt       = timestep;
  Src      = sourcestrength;
  rho      = density;
  mu       = viscosity;
  adhesionmagnitude = adhesion;
  ctension = tension;
  r_e      = re;
  xg       = cube(Ng+1,Ng+1,2,fill::zeros);
  sg       = mat(Ng+1,Ng+1,fill::zeros);
  fg       = cube(Ng+1,Ng+1,2,fill::zeros);
  vg       = cube(Ng+1,Ng+1,2,fill::zeros);
  ug       = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal = mat(2,0,fill::zeros);
  ubglobal = mat(2,0,fill::zeros);
  fbglobal = mat(2,0,fill::zeros);
  Nbs      = 4;
  xmin     =cen-static_cast<float>(dimensions);
  xmax     =cen+static_cast<float>(dimensions);
  sb       = mat(2,Nbs,fill::zeros);
  sbb      = mat(2,Nbs,fill::zeros);
  sb(0,0)  = 0.0; // One sink point at the centre of each fluid grid edge
  sb(1,0)  = xmin;
  sb(0,1)  = 0.0;
  sb(1,1)  = xmax;
  sb(0,2)  = xmin;
  sb(1,2)  = 0.0;
  sb(0,3)  = xmax;
  sb(1,3)  = 0.0;
  vector<element> Elements;
  h = cellheight;

  hg       = float(xmax-xmin)/float(Ng);
  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg+hg/2;
      xg(ii,jj,1)=xmin+jj*hg+hg/2;
    }
  }
}

void tissue::UpdatePositions(){
  xbglobal = xbglobal + dt*ubglobal;
}

tissue::~tissue() {}
