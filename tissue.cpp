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

tissue::tissue(const int& GridSize,const int& dimensions,const int& boundarypoints,const float& sourcestrength){
  Ng       = GridSize;
  Nbcell   = boundarypoints;
  Nb       = 0;
  Src      = sourcestrength;
  xg       = cube(Ng+1,Ng+1,2,fill::zeros);
  sg       = mat(Ng+1,Ng+1,fill::zeros);
  fg       = cube(Ng+1,Ng+1,2,fill::zeros);
  vg       = cube(Ng+1,Ng+1,2,fill::zeros);
  ug       = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal = mat(2,0,fill::zeros);
  ubglobal = mat(2,Nb,fill::zeros);
  fbglobal = mat(2,Nb,fill::zeros);
//  vector<cell> Cells;
  Nc  = 0;
  Nbs = 4;
  xmin=-static_cast<float>(dimensions);
  xmax=static_cast<float>(dimensions);
  sb  = mat(2,Nbs,fill::zeros);
  sbb = mat(2,Nbs,fill::zeros);
  // One sink point at the centre of each fluid grid edge
  sb(0,0) = 0.0;
  sb(1,0) = xmin;
  sb(0,1) = 0.0;
  sb(1,1) = xmax;
  sb(0,2) = xmin;
  sb(1,2) = 0.0;
  sb(0,3) = xmax;
  sb(1,3) = 0.0;

  hg  =float(xmax-xmin)/float(Ng);
  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg;
      xg(ii,jj,1)=xmin+jj*hg;
    }
  }
}

void tissue::CombineBoundaries(void){
  if (xbglobal.n_cols < Nc*Nbcell){
    xbglobal.resize(2,Nc*Nbcell);
    fbglobal.resize(2,Nc*Nbcell);
    ubglobal.resize(2,Nc*Nbcell);
  }
  for (int ii=0;ii<Nc;ii++){
    xbglobal.cols(ii*Nbcell,(ii+1)*Nbcell-1) = Cells[ii].xb;
    fbglobal.cols(ii*Nbcell,(ii+1)*Nbcell-1) = Cells[ii].fb;
  }
}

void tissue::UpdateSources(){
  // Resize arrays if more sources have been added
  if(sb.n_cols<(Nbs)){
    sb.resize(2,Nbs);
    sbb.resize(2,Nbs);
  }
  for (int ii=0; ii<Nbs; ii++){
    if (ii<4){
      sbb(0,ii)= -Src*Nc/4;
    }else {
      sbb(0,ii)= Src;
      sb(0,ii) = Cells[ii-4].com(0);
      sb(1,ii) = Cells[ii-4].com(1);
    }
  }
}

void tissue::AddCell(const float& len, const float& initialx, const float& initialy) {
  Cells.push_back(cell(Nbcell,len,initialx,initialy));
  Nc++;
  Nb=Nb+Nbcell;
  Nbs++;
}

tissue::~tissue() {}
