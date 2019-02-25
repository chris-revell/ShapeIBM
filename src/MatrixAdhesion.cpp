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

  Cell.NormaliseAdhesion();

  for (int ii=0; ii<Cell.Nb; ii++){
    // Locally store references to element under consideration and its neighbours by evaluating pointers in Tissue.Elements array
    element& elementii = Cell.Elements[Cell.ElementLabels[ii]];
    elementii.SetAdhesion();
    R = elementii.pos-Cell.com;
    F = -elementii.fb;
    magF = sqrt(dot(F,F));
    FdotR = dot(F,R);
    if (magF>elementii.adhesionmagnitude){
      //F = F*elementii.adhesionmagnitude/magF;
      F(0) = 0;
      F(1) = 0;
    }//else if (FdotR<0){
      //F(0) = 0;
      //F(1) = 0;
    //}
    //cout << "normalisationfactor " << elementii.normalisationfactor << endl;
    //cout << "F " << endl;
    //F.print();
    //cout << "F*normalisationfactor " << endl;
    //(F*elementii.normalisationfactor).print();
    //cout << "elementii.fb 0" << endl;
    //elementii.fb.print();
    elementii.fb = elementii.fb + F;//*elementii.normalisationfactor;
    //cout << "elementii.fb 1" << endl;
    //elementii.fb.print();
  }
}
