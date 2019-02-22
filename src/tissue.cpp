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
//#include <random>
#include <element.hpp>

using namespace arma;
using namespace std;

tissue::tissue(const int& GridSize,const int& dimensions,const int& boundarypoints,const float& sourcestrength,const float& density,const float& viscosity,const float& timestep){
  Ng       = GridSize;
  Nbcell   = boundarypoints;
  Nb       = 0;
  dt       = timestep;
  Src      = sourcestrength;
  rho      = density;
  mu       = viscosity;
  xg       = cube(Ng+1,Ng+1,2,fill::zeros);
//xNb      = cube(Ng+1,Ng+1,Nbcell,fill::zeros);
  sg       = mat(Ng+1,Ng+1,fill::zeros);
  fg       = cube(Ng+1,Ng+1,2,fill::zeros);
  vg       = cube(Ng+1,Ng+1,2,fill::zeros);
  ug       = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal = mat(2,0,fill::zeros);
  ubglobal = mat(2,0,fill::zeros);
  fbglobal = mat(2,0,fill::zeros);
  indices  = mat(2,0,fill::zeros);
  Nc       = 0;
  Nbs      = 4;
  xmin     =-static_cast<float>(dimensions);
  xmax     =static_cast<float>(dimensions);
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
  vector<cell> Cells;

  hg       = float(xmax-xmin)/float(Ng);
  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg+hg/2;
      xg(ii,jj,1)=xmin+jj*hg+hg/2;
    }
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

void tissue::AddCell(const float& len, const float& initialx, const float& initialy,const float& tension,const float& adhesion){
  Cells.push_back(cell(Nc,Nb,Nbcell,len,initialx,initialy,hg,tension,adhesion));
  Nc++;
  Nb=Nb+Nbcell;
  Nbs++;
}

//void tissue::BoundaryRefinement(){
//  float dx1,dy1,r1,newposx,newposy;
//  for (int kk=0;kk<Nc;kk++){
//    for (int ii=0;ii<Cells[kk].Nb;ii++){
//      dx1 = Cells[kk].Elements[Cells[kk].Elements[ii].neighbours[1]].pos(0)-Cells[kk].Elements[ii].pos(0);
//      dy1 = Cells[kk].Elements[Cells[kk].Elements[ii].neighbours[1]].pos(1)-Cells[kk].Elements[ii].pos(1);
//      // Find separation distances from x and y values.
//      r1=sqrt(pow(dx1,2)+pow(dy1,2));
//      if (r1>hg){
//        newposx = Cells[kk].Elements[ii].pos(0)+0.5*dx1;
//        newposy = Cells[kk].Elements[ii].pos(1)+0.5*dy1;
//        Cells[kk].Elements.push_back(element(Cells[kk].Elements[ii].ub(0),Cells[kk].Elements[ii].ub(1),kk,Nb,newposx,newposy,ii,((Cells[kk].Elements[ii].neighbours[1])%(Cells[kk].Nb+1)+(Cells[kk].Nb+1))%(Cells[kk].Nb+1)));
//        Cells[kk].Elements[ii].neighbours[1]=Cells[kk].Nb;
//        Nb++;
//        Cells[kk].Nb = Cells[kk].Nb+1;
//      }else if (r1<0.5*hg){
//
//      }
//    }
//  }
//}

void tissue::UpdatePositions(){
  xbglobal = xbglobal + dt*ubglobal;
}

tissue::~tissue() {}
