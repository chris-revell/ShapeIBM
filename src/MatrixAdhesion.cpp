//
//  MatrixAdhesion.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#include "MatrixAdhesion.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include <math.h>
#include "element.hpp"
#include <vector>

using namespace std;
using namespace arma;

void MatrixAdhesion(tissue& Tissue){

  //int xindex,yindex;
  float dr,Fmag,r0,r1,normalisationfactor;
  vec dvec = vec(2,fill::zeros);
  vec d0 = vec(2,fill::zeros);
  vec d1 = vec(2,fill::zeros);

  for (int ii=0; ii<Tissue.Nb; ii++){
    element& elementii = Tissue.Elements[ii];
    element& n0 = Tissue.Elements[elementii.n0];
    element& n1 = Tissue.Elements[elementii.n1];

    d0 = elementii.pos-n0.pos;
    d1 = elementii.pos-n1.pos;
    r0=sqrt(dot(d0,d0));
    r1=sqrt(dot(d1,d1));
    normalisationfactor = (r0+r1)/2.0;
    //xindex = floor((Tissue.Elements[ii].pos(0)-Tissue.xmin)/Tissue.hg);
    //yindex = floor((Tissue.Elements[ii].pos(1)-Tissue.xmin)/Tissue.hg);
    //dvec = vec(Tissue.xg(span(xindex),span(yindex),span::all))-Tissue.Elements[ii].pos;
    dvec = Tissue.Elements[ii].initialpos-Tissue.Elements[ii].pos;
    dr = sqrt(dot(dvec,dvec));
    if (dr > 0.000000001 && elementii.pos(0) < (2*Tissue.hg-Tissue.h/2) && elementii.pos(0) > (-2*Tissue.hg-Tissue.h/2) && elementii.pos(1) < (10*Tissue.hg) && elementii.pos(1) > (-9*Tissue.hg)){
      Fmag = 10*normalisationfactor*Tissue.adhesionmagnitude/dr;
    }else if (dr < 0.000000001 || dr > 2*Tissue.hg){
      Fmag=0;
    }else{
      Fmag = normalisationfactor*Tissue.adhesionmagnitude/dr;
    }
    elementii.fb = elementii.fb + Fmag*dvec;
  }
}
