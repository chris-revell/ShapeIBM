//
//  BoundToGrid1.cpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//

#include "BoundToGrid1.hpp"
#include <armadillo>
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

void BoundToGrid1(mat& sg,const int& Nbs,const mat& sb,const mat& sbb,const float& xmin,const float& xmax,const float& hg,const int& Ng){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  int Nx,Ny;

  sg.zeros();

  for (int n3=0;n3<Nbs;n3++){
    // Move points into the domain. xbb0 and xbb1 are the coordinates of the point in the square domain.
    xbb0=sb(0,n3);                                                 //IntoDom(sb(0,n3),xmin,xmax);
    xbb1=sb(1,n3);                                                 //IntoDom(sb(1,n3),xmin,xmax);
    // determine indices of the nearest lower-down grid point
    Nx=floor((xbb0-xmin)/hg);
    Ny=floor((xbb1-xmin)/hg);
    // tests all 16 possible grid points
    for (int ii=-1; ii<3; ii++){
      for (int jj=-1; jj<3; jj++){
        // compute the interpolation Delta function
        llx=xmin+Nx*hg+ii*hg;
        rr=abs(xbb0-llx);
        dx=DeltaFun(rr,0.5*hg);
        lly=xmin+Ny*hg+jj*hg;
        rr=abs(xbb1-lly);
        dy=DeltaFun(rr,0.5*hg);
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);
        // update the values if points are not passive
        if (dx*dy > 0){
          sg(x1,y1)  = sg(x1,y1) + sbb(0,n3)*dx*dy*0.5*hg;
          if (x2 != pas){
            sg(x2,y1)= sg(x2,y1) + sbb(0,n3)*dx*dy*0.5*hg;
          }
          if (y2 != pas){
            sg(x1,y2)= sg(x1,y2) + sbb(0,n3)*dx*dy*0.5*hg;
          }
          if ((x2 != pas) & (y2 != pas)){
            sg(x2,y2)= sg(x2,y2) + sbb(0,n3)*dx*dy*0.5*hg;
          }
        }
      }  // for jj
    } // for ii
  } // for n3
} // function BoundToGrid1
//-------------------------------------------------------------------//
