//
//  NavierStokes.cpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//

#include "NavierStokes.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "cell.hpp"
#include "smallfunctions.hpp"

using namespace arma;
using namespace std;


//---------------------------------------------------------------------//
// Uses the Fast Fourier Method to solve the Navier-Stokes equations   //
// for updating grid velocity ug due to grid forces fg and grid source //
// distribution sg, rho and mu are fluid constants, dt is a time step, //
// hg is a mesh width.                                                 //
//---------------------------------------------------------------------//
void NavierStokes(cube& vg,const cube& ug,const cube& fg,const mat& sg,const int& Ng,const float& rho,const float& mu,const float& dt,const float& hg){

    float pom,B1,B2,Aa,Bb,Bv;
    int in1,in2;
    float Eps=0.0000001;
    mat dummymat;
    cx_mat fvgg     = cx_mat(Ng,Ng,fill::zeros);
    cx_mat fsg      = cx_mat(Ng,Ng,fill::zeros);
    cx_mat fug1     = cx_mat(Ng,Ng,fill::zeros);
    cx_mat fug2     = cx_mat(Ng,Ng,fill::zeros);
    cube fvg0     = cube(Ng,Ng,2,fill::zeros);
    cube fvg1     = cube(Ng,Ng,2,fill::zeros);
    cx_mat vg0      = cx_mat(Ng,Ng,fill::zeros);
    cx_mat vg1      = cx_mat(Ng,Ng,fill::zeros);

    // stage n terms: force density fg, source distribution sg and current
    // velocity ug
    for (int n1=0; n1<Ng+1; n1++){
        for (int n2=0; n2<Ng+1; n2++){
            for (int ik=0; ik<2; ik++){
                // upwind scheme for the advection term
                if (ug(n1,n2,0) < 0){
                    in1=PeriodInd(n1,Ng,1);
                    pom=ug(in1,n2,ik)-ug(n1,n2,ik);
                }else{
                    in1=PeriodInd(n1,Ng,-1);
                    pom=ug(n1,n2,ik)-ug(in1,n2,ik);
                }
                vg(n1,n2,ik)=ug(n1,n2,0)*pom;

                if (ug(n1,n2,1) < 0){
                    in2=PeriodInd(n2,Ng,1);
                    pom=ug(n1,in2,ik)-ug(n1,n2,ik);
                }else{
                    in2=PeriodInd(n2,Ng,-1);
                    pom=ug(n1,n2,ik)-ug(n1,in2,ik);
                }
                vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,1)*pom;
                vg(n1,n2,ik)=-dt*vg(n1,n2,ik)/hg;

                // central difference for the grad of source term
                if (ik == 0){
                    in1=PeriodInd(n1,Ng,1);
                    in2=PeriodInd(n1,Ng,-1);
                    pom=sg(in1,n2)-sg(in2,n2);
                }else if (ik == 1){
                    in1=PeriodInd(n2,Ng,1);
                    in2=PeriodInd(n2,Ng,-1);
                    pom=sg(n1,in1)-sg(n1,in2);
                }
                vg(n1,n2,ik)=vg(n1,n2,ik)+dt*mu*pom/(6*hg*rho*rho);
                // current vlocity and force terms
                vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,ik)+dt*fg(n1,n2,ik)/rho;
            } // for ik
        } // for n2
    } // for n1

    // the Fast Fourier transforms of source distribution sg and stage n term vg
    fsg  = fft2(sg(span(0,Ng-1),span(0,Ng-1)));
    dummymat = vg.slice(0);
    fug1 = fft2(dummymat(span(0,Ng-1),span(0,Ng-1)));
    dummymat = vg.slice(1);
    fug2 = fft2(dummymat(span(0,Ng-1),span(0,Ng-1)));
    // determines fug - the Fourier Transform of the velocity field at the stage n+1
    for (int n1=0;n1<Ng;n1++){
        for (int n2=0;n2<Ng;n2++){

            B1=sin(2*M_PI*n1/Ng);
            B2=sin(2*M_PI*n2/Ng);
            Bb=pow(B1,2)+pow(B2,2);
            Aa=1+4*mu*dt*(pow(sin(M_PI*n1/Ng),2)+pow(sin(M_PI*n2/Ng),2))/(rho*hg*hg);

            if (Bb < Eps){
                fvg0(n1,n2,0)=real(fug1(n1,n2))/Aa;
                fvg0(n1,n2,1)=imag(fug1(n1,n2))/Aa;
                fvg1(n1,n2,0)=real(fug2(n1,n2))/Aa;
                fvg1(n1,n2,1)=imag(fug2(n1,n2))/Aa;
            }else{
                Bv=B1*real(fug1(n1,n2))+B2*real(fug2(n1,n2));
                fvg0(n1,n2,0)=(Bb*real(fug1(n1,n2))-B1*Bv)/(Aa*Bb);
                fvg1(n1,n2,0)=(Bb*real(fug2(n1,n2))-B2*Bv)/(Aa*Bb);

                Bv=B1*imag(fug1(n1,n2))+B2*imag(fug2(n1,n2));
                fvg0(n1,n2,1)=(Bb*imag(fug1(n1,n2))-B1*Bv)/(Aa*Bb);
                fvg1(n1,n2,1)=(Bb*imag(fug2(n1,n2))-B2*Bv)/(Aa*Bb);

                fvg0(n1,n2,0)=fvg0(n1,n2,0)+hg*B1*imag(fsg(n1,n2))/(Bb*rho);
                fvg0(n1,n2,1)=fvg0(n1,n2,1)-hg*B1*real(fsg(n1,n2))/(Bb*rho);
                fvg1(n1,n2,0)=fvg1(n1,n2,0)+hg*B2*imag(fsg(n1,n2))/(Bb*rho);
                fvg1(n1,n2,1)=fvg1(n1,n2,1)-hg*B2*real(fsg(n1,n2))/(Bb*rho);
            }
        } // for n2
    } // for n1

    // the inverse Fast Fourier Method of fvg
    dummymat = fvg0.slice(0);
    fvgg.set_real(dummymat(span(0,Ng-1),span(0,Ng-1)));
    dummymat = fvg0.slice(1);
    fvgg.set_imag(dummymat(span(0,Ng-1),span(0,Ng-1)));
    vg0 = ifft2(fvgg);

    dummymat = fvg1.slice(0);
    fvgg.set_real(dummymat(span(0,Ng-1),span(0,Ng-1)));
    dummymat = fvg1.slice(1);
    fvgg.set_imag(dummymat(span(0,Ng-1),span(0,Ng-1)));
    vg1 = ifft2(fvgg);

    for (int ii=0; ii<Ng; ii++){
        for (int jj=0; jj<Ng; jj++){
            vg(ii,jj,0)=real(vg0(ii,jj));
            vg(ii,jj,1)=real(vg1(ii,jj));
        }
    }

    for (int ii=0; ii<Ng; ii++){
        vg(Ng,ii,0)=vg(0,ii,0);
        vg(ii,Ng,0)=vg(ii,Ng,0);
        vg(Ng,ii,1)=vg(0,ii,1);
        vg(ii,Ng,1)=vg(ii,Ng,1);
    }
} // function NavierStokes
//-----------------------------------------------------------------------//
