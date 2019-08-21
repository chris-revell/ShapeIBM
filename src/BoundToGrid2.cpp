//
//  BoundToGrid2.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#include "BoundToGrid2.hpp"
#include <armadillo>
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

//-------------------------------------------------------------------//
// spreads the material values sb(1,Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1,2) in the //
// square domain (xmin,xmax)^2 with mesh width hg, material points   //
// separation hb and a radius of the discrete delta function hdl.    //
//-------------------------------------------------------------------//
void BoundToGrid2(cube& fg,const int& Nb,const mat& fbglobal,const mat& xbglobal,const float& xmin,const float& xmax,const float& hg,const int& Ng){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  int Nx,Ny;

  fg.zeros();

  for (int n3=0; n3<Nb; n3++){
    // move points into the domain
    xbb0=xbglobal(0,n3);                                     //IntoDom(xbglobal(0,n3),xmin,xmax);
    xbb1=xbglobal(1,n3);                                     //IntoDom(xbglobal(1,n3),xmin,xmax);
    // determine indices of the nearest lower-down grid point
    Nx=floor((xbb0-xmin)/hg);
    Ny=floor((xbb1-xmin)/hg);
    // tests all 16 possible grid points
    for (int ii=-1; ii<3;ii++){
      for (int jj=-1;jj<3;jj++){
        // compute the interpolation Delta function
        llx=xmin+Nx*hg+ii*hg;
        rr= abs(xbb0-llx);
        dx= DeltaFun(rr,hg);
        lly=xmin+Ny*hg+jj*hg;
        rr= abs(xbb1-lly);
        dy= DeltaFun(rr,hg);
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);
        // update the values if points are not pasive
        if (dx*dy > 0){
          fg(x1,y1,0)=fg(x1,y1,0)+fbglobal(0,n3)*dx*dy*0.5*hg;
          fg(x1,y1,1)=fg(x1,y1,1)+fbglobal(1,n3)*dx*dy*0.5*hg;
          if (x2 != pas){
            fg(x2,y1,0)=fg(x2,y1,0)+fbglobal(0,n3)*dx*dy*0.5*hg;
            fg(x2,y1,1)=fg(x2,y1,1)+fbglobal(1,n3)*dx*dy*0.5*hg;
          }
          if (y2 != pas){
            fg(x1,y2,0)=fg(x1,y2,0)+fbglobal(0,n3)*dx*dy*0.5*hg;
            fg(x1,y2,1)=fg(x1,y2,1)+fbglobal(1,n3)*dx*dy*0.5*hg;
          }
          if ((x2 != pas) && (y2 != pas)){
            fg(x2,y2,0)=fg(x2,y2,0)+fbglobal(0,n3)*dx*dy*0.5*hg;
            fg(x2,y2,1)=fg(x2,y2,1)+fbglobal(1,n3)*dx*dy*0.5*hg;
          }
        }
      } // jj
    } // ii
  } // for n3
} // function BoundToGrid2
//---------------------------------------------------------------------//
