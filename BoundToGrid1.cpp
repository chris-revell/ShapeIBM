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

//-----------------------------------------------------------------//
// spreads the material values sb(Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1) in the //
// square domain (xmin,xmax)^2 with mesh width hg, material points //
// separation hb and a radius of the discrete delta function hdl.  //
//-----------------------------------------------------------------//
void BoundToGrid1(tissue Tissue,const mat& xb,const mat& sb,const int& Nb){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  float xbb0;
  float xbb1;
  int Nx,Ny;

  int Ng  = Tissue.Ng;
  float hdl = Tissue.hg;
  float hg  = Tissue.hg;
  float hb  = 0.5*Tissue.hg;
  float xmn = Tissue.xmin;
  float xmx = Tissue.xmax;
  mat sg = Tissue.sg;

  for (int n3=0;n3<Nb;n3++){
    // Move points into the domain. xbb0 and xbb1 are the coordinates of the point in the square domain.
    xbb0=IntoDom(xb(0,n3),xmn,xmx);
    xbb1=IntoDom(xb(1,n3),xmn,xmx);
    // determine indices of the nearest lower-down grid point
    Nx=1+floor((xbb0-xmn)/hg);
    Ny=1+floor((xbb1-xmn)/hg);
    // tests all 16 possible grid points
    for (int ii=-1; ii<3; ii++){
      for (int jj=-1; jj<3; jj++){
        // compute the interpolation Delta function
        llx=xmn+(Nx-1)*hg+ii*hg;
        rr=fabs(xbb0-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmn+(Ny-1)*hg+jj*hg;
        rr=fabs(xbb1-lly);
        dy=DeltaFun(rr,hdl);
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmn,xmx);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmn,xmx);
        // update the values if points are not passive
        if (dx*dy > 0){
          sg(x1,y1)  = sg(x1,y1) + sb(0,n3)*dx*dy*hb;
          if (x2 != pas){
            sg(x2,y1)= sg(x2,y1) + sb(0,n3)*dx*dy*hb;
          }
          if (y2 != pas){
            sg(x1,y2)= sg(x1,y2) + sb(0,n3)*dx*dy*hb;
          }
          if ((x2 != pas) & (y2 != pas)){
            sg(x2,y2)= sg(x2,y2) + sb(0,n3)*dx*dy*hb;
          }
        }
      }  // for jj
    } // for ii
  } // for n3
} // function BoundToGrid1
//-------------------------------------------------------------------//
