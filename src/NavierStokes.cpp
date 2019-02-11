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
// for updating grid velocity Tissue.ug due to grid forces Tissue.fg and grid source //
// distribution Tissue.sg, rho and mu are fluid constants, dt is a time step, //
// hg is a mesh width.                                                 //
//---------------------------------------------------------------------//
void NavierStokes(tissue& Tissue){
  float pom,B1,B2,Aa,Bb,Bv;
  int in1,in2;
  float Eps=0.0000001;
  mat dummymat;
  cx_mat fvgg     = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  cx_mat fsg      = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  cx_mat fug1     = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  cx_mat fug2     = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  cube fvg0       = cube(Tissue.Ng,Tissue.Ng,2,fill::zeros);
  cube fvg1       = cube(Tissue.Ng,Tissue.Ng,2,fill::zeros);
  cx_mat vg0      = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  cx_mat vg1      = cx_mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  //mat  stoch_xb;          // Array containing stochastic update values for element positions
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Variables implementing normal distribution for stochastic dissipation
  //default_random_engine generator(seed);         // Variables implementing normal distribution for stochastic dissipation
  //normal_distribution<double> distribution;// Variables implementing normal distribution for stochastic dissipation

  //arma_rng::set_seed_random();
  //stoch_xb = mat(Tissue.Ng,Tissue.Ng,fill::zeros);
  //distribution = normal_distribution<double>(0.0,1.0);


  // stage n terms: force density Tissue.fg, source distribution Tissue.sg and current
  // velocity Tissue.ug
  for (int n1=0; n1<Tissue.Ng+1; n1++){
    for (int n2=0; n2<Tissue.Ng+1; n2++){
      for (int ik=0; ik<2; ik++){
        // upwind scheme for the advection term
        if (Tissue.ug(n1,n2,0) < 0){
          in1=PeriodInd(n1,Tissue.Ng,1);
          pom=Tissue.ug(in1,n2,ik)-Tissue.ug(n1,n2,ik);
        }else{
          in1=PeriodInd(n1,Tissue.Ng,-1);
          pom=Tissue.ug(n1,n2,ik)-Tissue.ug(in1,n2,ik);
        }
        Tissue.vg(n1,n2,ik)=Tissue.ug(n1,n2,0)*pom;
        if (Tissue.ug(n1,n2,1) < 0){
          in2=PeriodInd(n2,Tissue.Ng,1);
          pom=Tissue.ug(n1,in2,ik)-Tissue.ug(n1,n2,ik);
        }else{
          in2=PeriodInd(n2,Tissue.Ng,-1);
          pom=Tissue.ug(n1,n2,ik)-Tissue.ug(n1,in2,ik);
        }
        Tissue.vg(n1,n2,ik)=Tissue.vg(n1,n2,ik)+Tissue.ug(n1,n2,1)*pom;
        Tissue.vg(n1,n2,ik)=-Tissue.dt*Tissue.vg(n1,n2,ik)/Tissue.hg;
        // central difference for the grad of source term
        if (ik == 0){
          in1=PeriodInd(n1,Tissue.Ng,1);
          in2=PeriodInd(n1,Tissue.Ng,-1);
          pom=Tissue.sg(in1,n2)-Tissue.sg(in2,n2);
        }else if (ik == 1){
          in1=PeriodInd(n2,Tissue.Ng,1);
          in2=PeriodInd(n2,Tissue.Ng,-1);
          pom=Tissue.sg(n1,in1)-Tissue.sg(n1,in2);
        }
        Tissue.vg(n1,n2,ik)=Tissue.vg(n1,n2,ik)+Tissue.dt*Tissue.mu*pom/(6*Tissue.hg*Tissue.rho*Tissue.rho);
        // current vlocity and force terms
        Tissue.vg(n1,n2,ik)=Tissue.vg(n1,n2,ik)+Tissue.ug(n1,n2,ik)+Tissue.dt*Tissue.fg(n1,n2,ik)/Tissue.rho;
      } // for ik
    } // for n2
  } // for n1

  // the Fast Fourier transforms of source distribution Tissue.sg and stage n term Tissue.vg
  fsg  = fft2(Tissue.sg(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  dummymat = Tissue.vg.slice(0);
  fug1 = fft2(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  dummymat = Tissue.vg.slice(1);
  fug2 = fft2(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  // determines fug - the Fourier Transform of the velocity field at the stage n+1
  for (int n1=0;n1<Tissue.Ng;n1++){
    for (int n2=0;n2<Tissue.Ng;n2++){
      B1=sin(2*M_PI*n1/Tissue.Ng);
      B2=sin(2*M_PI*n2/Tissue.Ng);
      Bb=pow(B1,2)+pow(B2,2);
      Aa=1+4*Tissue.mu*Tissue.dt*(pow(sin(M_PI*n1/Tissue.Ng),2)+pow(sin(M_PI*n2/Tissue.Ng),2))/(Tissue.rho*Tissue.hg*Tissue.hg);
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
        fvg0(n1,n2,0)=fvg0(n1,n2,0)+Tissue.hg*B1*imag(fsg(n1,n2))/(Bb*Tissue.rho);
        fvg0(n1,n2,1)=fvg0(n1,n2,1)-Tissue.hg*B1*real(fsg(n1,n2))/(Bb*Tissue.rho);
        fvg1(n1,n2,0)=fvg1(n1,n2,0)+Tissue.hg*B2*imag(fsg(n1,n2))/(Bb*Tissue.rho);
        fvg1(n1,n2,1)=fvg1(n1,n2,1)-Tissue.hg*B2*real(fsg(n1,n2))/(Bb*Tissue.rho);
      }
    } // for n2
  } // for n1
  // the inverse Fast Fourier Method of fvg
  dummymat = fvg0.slice(0);
  fvgg.set_real(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  dummymat = fvg0.slice(1);
  fvgg.set_imag(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  vg0 = ifft2(fvgg);
  dummymat = fvg1.slice(0);
  fvgg.set_real(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  dummymat = fvg1.slice(1);
  fvgg.set_imag(dummymat(span(0,Tissue.Ng-1),span(0,Tissue.Ng-1)));
  vg1 = ifft2(fvgg);

  for (int ii=0; ii<Tissue.Ng; ii++){
    for (int jj=0; jj<Tissue.Ng; jj++){
      Tissue.vg(ii,jj,0)=real(vg0(ii,jj));
      Tissue.vg(ii,jj,1)=real(vg1(ii,jj));
    }
  }

  for (int ii=0; ii<Tissue.Ng; ii++){
    Tissue.vg(Tissue.Ng,ii,0)=Tissue.vg(0,ii,0);
    Tissue.vg(ii,Tissue.Ng,0)=Tissue.vg(ii,0,0);
    Tissue.vg(Tissue.Ng,ii,1)=Tissue.vg(0,ii,1);
    Tissue.vg(ii,Tissue.Ng,1)=Tissue.vg(ii,0,1);
  }
} 
