//
//  BoundToGrid1.cpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//

#include "BoundToGrid1.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "cell.hpp"
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

void BoundToGrid1(tissue& Tissue){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  int Nx,Ny;

  Tissue.sg.zeros();

  for (int n3=0;n3<Tissue.Nbs;n3++){
    // Move points into the domain. xbb0 and xbb1 are the coordinates of the point in the square domain.
    xbb0=Tissue.sb(0,n3);                                                 //IntoDom(Tissue.sb(0,n3),Tissue.xmin,Tissue.xmax);
    xbb1=Tissue.sb(1,n3);                                                 //IntoDom(Tissue.sb(1,n3),Tissue.xmin,Tissue.xmax);
    // determine indices of the nearest lower-down grid point
    Nx=floor((xbb0-Tissue.xmin)/Tissue.hg);
    Ny=floor((xbb1-Tissue.xmin)/Tissue.hg);
    // tests all 16 possible grid points
    for (int ii=-1; ii<3; ii++){
      for (int jj=-1; jj<3; jj++){
        // compute the interpolation Delta function
        llx=Tissue.xmin+Nx*Tissue.hg+ii*Tissue.hg;
        rr=abs(xbb0-llx);
        dx=DeltaFun(rr,0.5*Tissue.hg);
        lly=Tissue.xmin+Ny*Tissue.hg+jj*Tissue.hg;
        rr=abs(xbb1-lly);
        dy=DeltaFun(rr,0.5*Tissue.hg);
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        IndDel(y1,y2,lly,jj,Ny,Tissue.Ng,Tissue.xmin,Tissue.xmax);
        // update the values if points are not passive
        if (dx*dy > 0){
          Tissue.sg(x1,y1)  = Tissue.sg(x1,y1) + Tissue.sbb(0,n3)*dx*dy*0.5*Tissue.hg;
          if (x2 != pas){
            Tissue.sg(x2,y1)= Tissue.sg(x2,y1) + Tissue.sbb(0,n3)*dx*dy*0.5*Tissue.hg;
          }
          if (y2 != pas){
            Tissue.sg(x1,y2)= Tissue.sg(x1,y2) + Tissue.sbb(0,n3)*dx*dy*0.5*Tissue.hg;
          }
          if ((x2 != pas) & (y2 != pas)){
            Tissue.sg(x2,y2)= Tissue.sg(x2,y2) + Tissue.sbb(0,n3)*dx*dy*0.5*Tissue.hg;
          }
        }
      }  // for jj
    } // for ii
  } // for n3
} // function BoundToGrid1
//-------------------------------------------------------------------//
