//
//  GridToBound.cpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//

#include "GridToBound.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

//-----------------------------------------------------------------------//
// interpolates the grid values Tissue.vg(Tissue.Ng,Tissue.Ng,2) (velocities) to the material //
// values Tissue.ubglobal(1,Nb) defined at the material points xb(1,Nb) in the square //
// domain (xmn,xmx)^2 with mesh width hg and a radius of the discrete    //
// Dirac delta Tissue.hg.                                                      //
//-----------------------------------------------------------------------//
void GridToBound(tissue& Tissue){

  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int Nx,Ny;
  int x1 = 0;
  int x2 = 0;
  int y1 = 0;
  int y2 = 0;
  for (int n3=0; n3<Tissue.Nb; n3++){
    // moves points into the domain
    xbb0=Tissue.xbglobal(0,n3);                                                   //IntoDom(Tissue.xbglobal(0,n3),Tissue.xmin,Tissue.xmax);
    xbb1=Tissue.xbglobal(1,n3);                                                   //IntoDom(Tissue.xbglobal(1,n3),Tissue.xmin,Tissue.xmax);
    Nx=floor((xbb0-Tissue.xmin)/Tissue.hg); // Indices of grid points containing boundary point n3
    Ny=floor((xbb1-Tissue.xmin)/Tissue.hg); // Indices of grid points containing boundary point n3
    // test all 16 neighboring points
    for (int ii=-1;ii<3;ii++){
      for (int jj=-1;jj<3;jj++){
        // determine the value of the interpolation Delta-function
        llx=Tissue.xmin+Nx*Tissue.hg+ii*Tissue.hg; // x dimension coordinate of lower x edge of neighbouring fluid grid point
        rr=abs(xbb0-llx);// x dimension distance from boundary element n3 to the lower x edge of the neighbouring fluid grid point
        dx=DeltaFun(rr,Tissue.hg); // Value of spread function at this x distance from boundary point n3
        lly=Tissue.xmin+Ny*Tissue.hg+jj*Tissue.hg; // y dimension coordinate of lower y edge of neighbouring fluid grid point
        rr=abs(xbb1-lly);// y dimension distance from boundary element n3 to the lower y edge of the neighbouring fluid grid point
        dy=DeltaFun(rr,Tissue.hg); // Value of spread function at this x distance from boundary point n3
        // determine the indices of grid points gaining positive impact
        IndDel(x1,x2,llx,ii,Nx,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        IndDel(y1,y2,lly,jj,Ny,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        // update the values if inside the impact domain
        if (dx*dy > 0){
          Tissue.ubglobal(0,n3)=Tissue.ubglobal(0,n3)+Tissue.vg(x1,y1,0)*dx*dy*Tissue.hg*Tissue.hg;
          Tissue.ubglobal(1,n3)=Tissue.ubglobal(1,n3)+Tissue.vg(x1,y1,1)*dx*dy*Tissue.hg*Tissue.hg;
        }
      } // for jj
    } // for ii
  }// for n3
} // function GridToBound
//-------------------------------------------------------------------//
