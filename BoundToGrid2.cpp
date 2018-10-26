//
//  BoundToGrid2.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#include "BoundToGrid2.hpp"

//-------------------------------------------------------------------//
// spreads the material values sb(1,Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1,2) in the //
// square domain (xmin,xmax)^2 with mesh width hg, material points   //
// separation hb and a radius of the discrete delta function hdl.    //
//-------------------------------------------------------------------//
void BoundToGrid2(cube& sg,const mat& xb,const mat& sb,const int& Nb,const int& Ng,const float& hdl,const float& hg,const float& hb,const float& xmin,const float& xmax){

    // passive value - do nothing
    int pas=-100;
    int Nx,Ny;
    float xbb0,xbb1;
    int x1=0;
    int x2=0;
    int y1=0;
    int y2=0;
    float llx,rr,dx,lly,dy;

    for (int n3=0; n3<Nb; n3++){
        // move points into the domain

        xbb0=IntoDom(xb(0,n3),xmin,xmax);

        xbb1=IntoDom(xb(1,n3),xmin,xmax);

        // determine indices of the nearest lower-down grid point
        Nx=1+floor((xbb0-xmin)/hg);
        Ny=1+floor((xbb1-xmin)/hg);

        // tests all 16 possible grid points
        for (int ii=-1; ii<2;ii++){
            for (int jj=-1;jj<2;jj++){
                // compute the interpolation Delta function
                llx=xmin+(Nx-1)*hg+ii*hg;
                rr= fabs(xbb0-llx);
                dx= DeltaFun(rr,hdl);
                lly=xmin+(Ny-1)*hg+jj*hg;
                rr= fabs(xbb1-lly);
                dy= DeltaFun(rr,hdl);

                // determine indices of the grid points to update
                IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
                IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);

                // update the values if points are not pasive
                if (dx*dy > 0){
                    sg(x1,y1,0)=sg(x1,y1,0)+sb(0,n3)*dx*dy*hb;
                    sg(x1,y1,1)=sg(x1,y1,1)+sb(1,n3)*dx*dy*hb;

                    if (x2 != pas){
                        sg(x2,y1,0)=sg(x2,y1,0)+sb(0,n3)*dx*dy*hb;
                        sg(x2,y1,1)=sg(x2,y1,1)+sb(1,n3)*dx*dy*hb;
                    }
                    if (y2 != pas){
                        sg(x1,y2,0)=sg(x1,y2,0)+sb(0,n3)*dx*dy*hb;
                        sg(x1,y2,1)=sg(x1,y2,1)+sb(1,n3)*dx*dy*hb;
                    }
                    if ((x2 != pas) && (y2 != pas)){
                        sg(x2,y2,0)=sg(x2,y2,0)+sb(0,n3)*dx*dy*hb;
                        sg(x2,y2,1)=sg(x2,y2,1)+sb(1,n3)*dx*dy*hb;
                    }
                }

            } // jj
        } // ii
    } // for n3
} // function BoundToGrid2
//---------------------------------------------------------------------//
