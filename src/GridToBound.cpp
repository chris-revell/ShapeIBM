//
//  GridToBound.cpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//

#include "GridToBound.hpp"
#include <armadillo>
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

//-----------------------------------------------------------------------//
// interpolates the grid values vg(Ng,Ng,2) (velocities) to the material //
// values ubglobal(1,Nb) defined at the material points xb(1,Nb) in the square //
// domain (xmn,xmx)^2 with mesh width hg and a radius of the discrete    //
// Dirac delta hg.                                                      //
//-----------------------------------------------------------------------//
void GridToBound(mat& ubglobal,const mat& xbglobal,const cube& vg,const float& hg,const float& xmin,const float& xmax,const int& Nb,const int& Ng){

  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int Nx,Ny;
  int x1 = 0;
  int x2 = 0;
  int y1 = 0;
  int y2 = 0;
  for (int n3=0; n3<Nb; n3++){
    // moves points into the domain
    xbb0=xbglobal(0,n3);                                                   //IntoDom(xbglobal(0,n3),xmin,xmax);
    xbb1=xbglobal(1,n3);                                                   //IntoDom(xbglobal(1,n3),xmin,xmax);
    Nx=floor((xbb0-xmin)/hg); // Indices of grid points containing boundary point n3
    Ny=floor((xbb1-xmin)/hg); // Indices of grid points containing boundary point n3
    // test all 16 neighboring points
    for (int ii=-1;ii<3;ii++){
      for (int jj=-1;jj<3;jj++){
        // determine the value of the interpolation Delta-function
        llx=xmin+Nx*hg+ii*hg; // x dimension coordinate of lower x edge of neighbouring fluid grid point
        rr=abs(xbb0-llx);// x dimension distance from boundary element n3 to the lower x edge of the neighbouring fluid grid point
        dx=DeltaFun(rr,hg); // Value of spread function at this x distance from boundary point n3
        lly=xmin+Ny*hg+jj*hg; // y dimension coordinate of lower y edge of neighbouring fluid grid point
        rr=abs(xbb1-lly);// y dimension distance from boundary element n3 to the lower y edge of the neighbouring fluid grid point
        dy=DeltaFun(rr,hg); // Value of spread function at this x distance from boundary point n3
        // determine the indices of grid points gaining positive impact
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);
        // update the values if inside the impact domain
        if (dx*dy > 0){
          ubglobal(0,n3)=ubglobal(0,n3)+vg(x1,y1,0)*dx*dy*hg*hg;
          ubglobal(1,n3)=ubglobal(1,n3)+vg(x1,y1,1)*dx*dy*hg*hg;
        }
      } // for jj
    } // for ii
  }// for n3
} // function GridToBound
//-------------------------------------------------------------------//
