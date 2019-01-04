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
#include <random>

using namespace arma;
using namespace std;

tissue::tissue(const int& GridSize,const int& dimensions,const int& boundarypoints,const float& sourcestrength,const float& density,const float& viscosity,const float& stoch,const float& timestep){
  Ng       = GridSize;
  Nbcell   = boundarypoints;
  Nb       = 0;
  dt       = timestep;
  Src      = sourcestrength;
  rho      = density;
  mu       = viscosity;
  xi     = stoch;
  xg       = cube(Ng+1,Ng+1,2,fill::zeros);
  sg       = mat(Ng+1,Ng+1,fill::zeros);
  fg       = cube(Ng+1,Ng+1,2,fill::zeros);
  vg       = cube(Ng+1,Ng+1,2,fill::zeros);
  ug       = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal = mat(2,0,fill::zeros);
  ubglobal = mat(2,0,fill::zeros);
  fbglobal = mat(2,0,fill::zeros);
  indices  = mat(2,0,fill::zeros);
  stoch_xb = mat(2,0,arma::fill::zeros);
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
  distribution = std::normal_distribution<double>(0.0,2.0);

  hg  =float(xmax-xmin)/float(Ng);
  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg+hg/2;
      xg(ii,jj,1)=xmin+jj*hg+hg/2;
    }
  }
}

void tissue::CombineBoundaries(void){
  Nb = 0;
  for (int ii=0;ii<Nc;ii++){
    Nb = Nb+Cells[ii].Nb;
  }
  if (xbglobal.n_cols < Nb){
    xbglobal.resize(2,Nb);
    fbglobal.resize(2,Nb);
    ubglobal.resize(2,Nb);
    indices.resize(2,Nb);
    stoch_xb.resize(2,Nb);
  }
  for (int ii=0;ii<Nc;ii++){
    for (int jj=0; jj<Cells[ii].Nb;jj++){
      xbglobal(0,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].pos(0);
      xbglobal(1,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].pos(1);
      fbglobal(0,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].internalforce(0);
      fbglobal(1,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].internalforce(1);
      ubglobal(0,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].ub(0);
      ubglobal(1,Cells[ii].Elements[jj].label) = Cells[ii].Elements[jj].ub(1);
      indices(0,Cells[ii].Elements[jj].label) = ii;
      indices(1,Cells[ii].Elements[jj].label) = jj;
    }
  }
  stoch_xb.randn(); // Set random angle
  for (int ii=0;ii<Nb;ii++){
    float gauss = xi*distribution(generator);       // Find random magnitude
    stoch_xb(0,ii) = gauss*cos(2*M_PI*stoch_xb(0,ii));// Use random angle to find magnitude in each dimension
    stoch_xb(1,ii) = gauss*sin(2*M_PI*stoch_xb(1,ii));// Use random angle to find magnitude in each dimension
  }
  fbglobal = fbglobal+stoch_xb/1000.0;
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

void tissue::AddCell(const float& len, const float& initialx, const float& initialy,float& tension){
  Cells.push_back(cell(Nc,Nb,Nbcell,len,initialx,initialy,hg,tension));
  Nc++;
  Nb=Nb+Nbcell;
  Nbs++;
}

void tissue::BoundaryRefinement(){
  float dx1,dy1,r1,newposx,newposy;
  for (int kk=0;kk<Nc;kk++){
    for (int ii=0;ii<Cells[kk].Nb;ii++){

      dx1 = Cells[kk].Elements[Cells[kk].Elements[ii].neighbours[1]].pos(0)-Cells[kk].Elements[ii].pos(0);
      dy1 = Cells[kk].Elements[Cells[kk].Elements[ii].neighbours[1]].pos(1)-Cells[kk].Elements[ii].pos(1);
      // Find separation distances from x and y values.
      r1=sqrt(pow(dx1,2)+pow(dy1,2));
      if (r1>hg){
        newposx = Cells[kk].Elements[ii].pos(0)+0.5*dx1;
        newposy = Cells[kk].Elements[ii].pos(1)+0.5*dy1;
        Cells[kk].Elements.push_back(element(Cells[kk].Elements[ii].ub(0),Cells[kk].Elements[ii].ub(1),kk,Nb,newposx,newposy,ii,((Cells[kk].Elements[ii].neighbours[1])%(Cells[kk].Nb+1)+(Cells[kk].Nb+1))%(Cells[kk].Nb+1)));
        Cells[kk].Elements[ii].neighbours[1]=Cells[kk].Nb;
        Nb++;
        Cells[kk].Nb = Cells[kk].Nb+1;
      }
    }
  }
}

void tissue::UpdatePositions(){
  //stoch_xb.randn();
  xbglobal = xbglobal + dt*ubglobal;// + stoch_xb.tail_cols(xbglobal.n_cols-1)/1000;
  //xbglobal.tail_cols(xbglobal.n_cols-1) = xbglobal.tail_cols(xbglobal.n_cols-1) + dt*ubglobal.tail_cols(xbglobal.n_cols-1);// + stoch_xb.tail_cols(xbglobal.n_cols-1)/1000;
}

tissue::~tissue() {}
