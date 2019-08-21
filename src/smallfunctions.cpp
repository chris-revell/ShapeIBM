//
// Created by Christopher Revell on 24/10/2018.
// Declaration of IBM functions

#include "smallfunctions.hpp"
#include <math.h>

using namespace std;

//--------------------------------------------------------------------//
// transforms the real coordinate xy of the body into the correspon-  //
// ding coordinate pom inside the periodic domain (xmn,xmx)x(xmn,xmx) //                             //
//--------------------------------------------------------------------//
//float IntoDom(const float& xy,const float& xmin,const float& xmax){
//    float len;
//    float pom;
//    len=xmax-xmin;
//    pom=xy;
//    while (pom>xmax){
//        pom=pom-len;
//    }
//    while (pom<xmin){
//        pom=pom+len;
//    }
//    return pom;
//} // function IntoDom
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// determines value of the unitary bell-shaped discrete approximation //
// to the Dirac delta function of radius h for a point located at    //
// distance r from the center.                                       //
//-------------------------------------------------------------------//
float DeltaFun(const float& r,const float& h){
    float dist;

    if (abs(r) < (2*h)){
        dist=0.25*(1+cos(0.5*M_PI*r/h))/h;
    }else{
        dist=0;
    }
    return dist;
} // function DeltaFun
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// determines indices in1 and in2 of grid elements which need  to be //
// updated in the interpolation of point ll (shifted by ij according //
// to the reference point Nij). Depending on the location of point   //
// ll, there may be at most two active indeces for each element, due //
// to periodicity of fluid domain.                                   //
//   0 elements to update, if ll inside the domain,                  //
//   1 element  to update, if ll is close to the boundary,           //
//   2 elements to update, if ll is close to the corner.             //
//-------------------------------------------------------------------//
void IndDel(int& in1,int& in2,const float& ll,const int& ij,const int& Nij,const int& Nmx,const float& xmn,const float& xmx){
    // passive value - do nothing
    int pas=-100;

    in2=pas;
    if (ll < xmn){
        in1=Nmx+ij;
    }else if ((ll == xmn) || (ll == xmx)){
        in1=Nmx;
        in2=0;
    }else if (ll > xmx){
        in1=ij-1;
    }else{
        in1=Nij+ij-1;
    }
} // function IndDel
//--------------------------------------------------------------------//

//-------------------------------------------------------------------//
// computes shifting of the index nn in the square periodic domain   //
// of size Ng+1 where the 1st and the (Ng+1)st columns and rows are  //
// identical according to the value of which:                        //
//    which = 1 -- nn shifted one element to the right               //
//    which =-1 -- nn shifted one element to the left                //
//-------------------------------------------------------------------//
int PeriodInd(const int& nn,const int& Ng,const int& which){
    int ind=0;
    if (which == (-1)){
        if (nn == 0){
            ind=Ng-1;
        }else{
            ind=nn-1;
        }
    }
    if (which == 1){
        if (nn == Ng){
            ind=1;
        }else{
            ind=nn+1;
        }
    }
    return ind;
} // function PeriodInd
//------------------------------------------------------------------//
