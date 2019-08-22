//
//  MatrixAdhesion.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#include "MatrixAdhesion.hpp"
#include <armadillo>
#include <vector>
#include <math.h>
#include "element.hpp"

using namespace std;
using namespace arma;

void MatrixAdhesion(vector<element>& Elements,const int& Nb,const float& adhesion,const float& hg){

  //int xindex,yindex;
  float dr,Fmag,r0,r1,normalisationfactor;
  vec dvec = vec(2,fill::zeros);
  vec d0 = vec(2,fill::zeros);
  vec d1 = vec(2,fill::zeros);

  for (int ii=0; ii<Nb; ii++){
    element& elementii = Elements[ii];
    element& n0 = Elements[elementii.n0];
    element& n1 = Elements[elementii.n1];

    d0 = elementii.pos-n0.pos;
    d1 = elementii.pos-n1.pos;
    r0=sqrt(dot(d0,d0));
    r1=sqrt(dot(d1,d1));
    normalisationfactor = (r0+r1)/2.0;
    dvec = Elements[ii].initialpos-Elements[ii].pos;
    dr = sqrt(dot(dvec,dvec));
    if (dr < 0.000000001 || dr > 2*hg){
      Fmag=0;
    }else{
      Fmag = normalisationfactor*adhesion;
    }
    elementii.fb = elementii.fb + Fmag*dvec;
  }
}
