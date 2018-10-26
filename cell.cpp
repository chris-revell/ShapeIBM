//
// Created by Christopher Revell on 24/10/2018.
//

#include "cell.hpp"
#include "armadillo"
#include "cmath"

cell::cell(const int& Nb, const float& len, const float& initialx, const float& initialy) { // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
    hb=2*M_PI/static_cast<float>(Nb); // Typical angular spacing between boundary elements given typical radius len
    xb = arma::mat(2,Nb,arma::fill::zeros); // Positions of all boundary points in cell
    fb = arma::mat(2,Nb,arma::fill::zeros); // Forces on all boundary points in cell arising from interactions with other boundary points
    corticaltension = 1; // Spring constant of boundary forces
    // Loop over initial element angles to set cartesian coordinates of all boundary points. Add constant value to set initial cell position.
    for (int ii=0;ii<Nb;ii++){
      xb(0,ii) = len*cos(ii*hb)+initialx;
      xb(1,ii) = len*sin(ii*hb)+initialy;
    }
}

void AdjacentForces(); {
  float dl1,dl2,dr1,dr2,ndl,ndr;
  for (int ii=0; ii<Nb; ii++) {
      if (ii==0) {
          dl1=xb(0,Nb-1)-xb(0,ii);
          dl2=xb(1,Nb-1)-xb(1,ii);
          dr1=xb(0,ii+1)-xb(0,ii);
          dr2=xb(1,ii+1)-xb(1,ii);
      } else if (ii==Nb-1) {
          dl1=xb(0,ii-1)-xb(0,ii);
          dl2=xb(1,ii-1)-xb(1,ii);
          dr1=xb(0,0)-xb(0,ii);
          dr2=xb(1,0)-xb(1,ii);
      } else {
          dl1=xb(0,ii-1)-xb(0,ii);
          dl2=xb(1,ii-1)-xb(1,ii);
          dr1=xb(0,ii+1)-xb(0,ii);
          dr2=xb(1,ii+1)-xb(1,ii);
      }
      ndl=sqrt(pow(dl1,2)+pow(dl2,2));
      ndr=sqrt(pow(dr1,2)+pow(dr2,2));
      fb(0,ii)=corticaltension*((ndl-hb)*dl1/ndl+(ndr-hb)*dr1/ndr);
      fb(1,ii)=corticaltension*((ndl-hb)*dl2/ndl+(ndr-hb)*dr2/ndr);
  }
}



cell::~cell() {}
