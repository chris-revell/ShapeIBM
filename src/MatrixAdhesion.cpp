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
#include "element.hpp"
#include <math.h>

using namespace std;
using namespace arma;

void MatrixAdhesion(tissue& Tissue){

  int xindex,yindex;
  float dr,Fmag;
  vec dvec = vec(2,fill::zeros);

  for (int jj=0;jj<Tissue.Nc;jj++){
    cell& Cell = Tissue.Cells[jj];
    Cell.NormaliseAdhesion();
    for (int ii=0; ii<Cell.Nb; ii++){
      element& elementii = Cell.Elements[Cell.ElementLabels[ii]];
      elementii.SetAdhesion();
      xindex = floor((elementii.pos(0)-Tissue.xmin)/Tissue.hg);
      yindex = floor((elementii.pos(1)-Tissue.xmin)/Tissue.hg);
      dvec = vec(Tissue.xg(span(xindex),span(yindex),span::all))-elementii.pos;
      dr = sqrt(dot(dvec,dvec));
      //if (dr<Tissue.hg){
        Fmag = elementii.adhesionmagnitude/dr;//*dr/Tissue.hg;
        elementii.fb = elementii.fb + Fmag*dvec;
      //}else{
      //}
    }
  }
}
