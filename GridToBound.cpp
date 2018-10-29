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
#include "cell.hpp"
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;

//-----------------------------------------------------------------------//
// interpolates the grid values fg(Ng,Ng,2) (velocities) to the material //
// values fb(1,Nb) defined at the material points xb(1,Nb) in the square //
// domain (xmn,xmx)^2 with mesh width hg and a radius of the discrete    //
// Dirac delta hdl.                                                      //
//-----------------------------------------------------------------------//
void GridToBound(mat& fb,const mat& xb,const int& Nb,const cube& fg,const int& Ng,const float& hdl,const float& hg,const float& xmn,const float& xmx){

    float llx,rr,dx,lly,dy,xbb0,xbb1;
    int Nx,Ny;
    int x1 = 0;
    int x2 = 0;
    int y1 = 0;
    int y2 = 0;


    for (int n3=0; n3<Nb; n3++){
        // moves points into the domain
        xbb0=IntoDom(xb(0,n3),xmn,xmx);
        xbb1=IntoDom(xb(1,n3),xmn,xmx);
        //computes the indices of the nearest down-left grid point
        Nx=1+floor((xbb0-xmn)/hg);
        Ny=1+floor((xbb1-xmn)/hg);
        // test all 16 neighboring points
        for (int ii=-1;ii<3;ii++){
            for (int jj=-1;jj<3;jj++){
                // determine the value of the interpolation Delta-function
                llx=xmn+(Nx-1)*hg+ii*hg;
                rr=fabs(xbb0-llx);
                dx=DeltaFun(rr,hdl);
                lly=xmn+(Ny-1)*hg+jj*hg;
                rr=fabs(xbb1-lly);
                dy=DeltaFun(rr,hdl);

                // determine the indices of grid points gaining positive impact
                IndDel(x1,x2,llx,ii,Nx,Ng,xmn,xmx);
                IndDel(y1,y2,lly,jj,Ny,Ng,xmn,xmx);

                // update the values if inside the impact domain
                if (dx*dy > 0){
                    fb(0,n3)=fb(0,n3)+fg(x1,y1,0)*dx*dy*hg*hg;
                    fb(1,n3)=fb(1,n3)+fg(x1,y1,1)*dx*dy*hg*hg;
                }
            } // for jj
        } // for ii
    }// for n3
} // function GridToBound
//-------------------------------------------------------------------//
