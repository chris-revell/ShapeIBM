//
//  BoundToGrid2.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#include "BoundToGrid2.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "cell.hpp"
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

//-------------------------------------------------------------------//
// spreads the material values sb(1,Tissue.Nb) (forces, sources) defined at //
// material points xb(1,Tissue.Nb) to the fluid grid sg(Tissue.Ng+1,Tissue.Ng+1,2) in the //
// square domain (xmin,Tissue.xmax)^2 with mesh width hg, material points   //
// separation hb and a radius of the discrete delta function hdl.    //
//-------------------------------------------------------------------//
void BoundToGrid2(tissue& Tissue){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  int Nx,Ny;

  Tissue.fg.zeros();

  for (int n3=0; n3<Tissue.Nb; n3++){
    // move points into the domain
    xbb0=Tissue.xbglobal(0,n3);//IntoDom(Tissue.xbglobal(0,n3),Tissue.xmin,Tissue.xmax);
    xbb1=Tissue.xbglobal(1,n3);//IntoDom(Tissue.xbglobal(1,n3),Tissue.xmin,Tissue.xmax);
    // determine indices of the nearest lower-down grid point
    Nx=floor((xbb0-Tissue.xmin)/Tissue.hg);
    Ny=floor((xbb1-Tissue.xmin)/Tissue.hg);
    // tests all 16 possible grid points
    for (int ii=-1; ii<3;ii++){
      for (int jj=-1;jj<3;jj++){
        // compute the interpolation Delta function
        llx=Tissue.xmin+Nx*Tissue.hg+ii*Tissue.hg;
        rr= abs(xbb0-llx);
        dx= DeltaFun(rr,Tissue.hg);
        lly=Tissue.xmin+Ny*Tissue.hg+jj*Tissue.hg;
        rr= abs(xbb1-lly);
        dy= DeltaFun(rr,Tissue.hg);
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        IndDel(y1,y2,lly,jj,Ny,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        // update the values if points are not pasive
        if (dx*dy > 0){
          Tissue.fg(x1,y1,0)=Tissue.fg(x1,y1,0)+Tissue.fbglobal(0,n3)*dx*dy*0.5*Tissue.hg;
          Tissue.fg(x1,y1,1)=Tissue.fg(x1,y1,1)+Tissue.fbglobal(1,n3)*dx*dy*0.5*Tissue.hg;
          if (x2 != pas){
            Tissue.fg(x2,y1,0)=Tissue.fg(x2,y1,0)+Tissue.fbglobal(0,n3)*dx*dy*0.5*Tissue.hg;
            Tissue.fg(x2,y1,1)=Tissue.fg(x2,y1,1)+Tissue.fbglobal(1,n3)*dx*dy*0.5*Tissue.hg;
          }
          if (y2 != pas){
            Tissue.fg(x1,y2,0)=Tissue.fg(x1,y2,0)+Tissue.fbglobal(0,n3)*dx*dy*0.5*Tissue.hg;
            Tissue.fg(x1,y2,1)=Tissue.fg(x1,y2,1)+Tissue.fbglobal(1,n3)*dx*dy*0.5*Tissue.hg;
          }
          if ((x2 != pas) && (y2 != pas)){
            Tissue.fg(x2,y2,0)=Tissue.fg(x2,y2,0)+Tissue.fbglobal(0,n3)*dx*dy*0.5*Tissue.hg;
            Tissue.fg(x2,y2,1)=Tissue.fg(x2,y2,1)+Tissue.fbglobal(1,n3)*dx*dy*0.5*Tissue.hg;
          }
        }
      } // jj
    } // ii
  } // for n3
} // function BoundToGrid2
//---------------------------------------------------------------------//
