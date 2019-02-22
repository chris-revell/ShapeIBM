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

double SafeAcos (double x){
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acos (x) ;
  }
double safeatan2 (double x,double y){
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return atan2 (x,y) ;
  }

void MatrixAdhesion(cell& Cell){

  vec R = vec(2,fill::zeros);
  //vec A = vec(2,fill::zeros);
  //vec B = vec(2,fill::zeros);
  vec F = vec(2,fill::zeros);
  float magF,FdotR;//,theta,phi,gamma;

  for (int ii=0; ii<Cell.Nb; ii++){
    // Locally store references to element under consideration and its neighbours by evaluating pointers in Tissue.Elements array
    element& elementii = Cell.Elements[Cell.ElementLabels[ii]];
    //element& n0        = Cell.Elements[Cell.Elements[ii].neighbours[0]];
    //element& n1        = Cell.Elements[Cell.Elements[ii].neighbours[1]];
    // Find vectors from cell centre of mass to element and from neighbours to element.
    R = elementii.pos-Cell.com;
    //A = elementii.pos-n0.pos;
    //B = elementii.pos-n1.pos;
    //// Force is colinear with the sum of unit vectors from neighbours to element.
    //theta = M_PI - 0.5*SafeAcos(dot(A,B)/(sqrt(dot(A,A))*sqrt(dot(B,B))));
    //phi = safeatan2(B(0),B(1));
    //gamma = theta+phi+M_PI;
    //F(0) = sin(gamma);
    //F(1) = cos(gamma);
    //FdotR = dot(F,R);
    //if (FdotR<0){
    //  elementii.fb += ((-1.0)*F*elementii.adhesionmagnitude);
    //}else{
    //  elementii.fb += (F*elementii.adhesionmagnitude);
    //}
    F = -elementii.fb;
    magF = sqrt(dot(F,F));
    FdotR = dot(F,R);
    if (magF>elementii.adhesionmagnitude){
      F = F*elementii.adhesionmagnitude/magF;
    }else if (FdotR<0){
      F(0) = 0;
      F(1) = 0;
    }
    elementii.fb = elementii.fb + F;
  }
}
