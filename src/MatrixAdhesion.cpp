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
  vec A = vec(2,fill::zeros);
  vec B = vec(2,fill::zeros);
  vec F = vec(2,fill::zeros);
  float FdotR,magA,magB,magF,theta,phi,zeta,gamma;

  for (int ii=0; ii<Cell.Nb; ii++){
    // Locally store references to element under consideration and its neighbours by evaluating pointers in Tissue.Elements array
    element& elementii = Cell.Elements[ii];
    element& n0        = Cell.Elements[Cell.Elements[ii].neighbours[0]];
    element& n1        = Cell.Elements[Cell.Elements[ii].neighbours[1]];
    // Find vectors from cell centre of mass to element and from neighbours to element.
    R = elementii.pos-Cell.com;
    A = elementii.pos-n0.pos;
    B = elementii.pos-n1.pos;
    magA = sqrt(dot(A,A));
    magB = sqrt(dot(B,B));
    // Force is colinear with the sum of unit vectors from neighbours to element.
//    cout << "elementii.pos " << elementii.pos(0) << elementii.pos(1) << endl;
//    cout << "n0.pos " << n0.pos(0) << n0.pos(1) << endl;
//    cout << "n1.pos " << n1.pos(0) << n1.pos(1) << endl;
//    cout << "Cell.com " << Cell.com(0) << Cell.com(1) << endl;
//    cout << "A " << A(0) << " " << A(1) << endl;
//    cout << "B " << B(0) << " " << B(1) << endl;
//    cout << "R " << R(0) << " " << R(1) << endl;
//    cout << "dot(A,B) " << dot(A,B) << endl;
//    cout << "magA " << magA << endl;
//    cout << "magB " << magB << endl;
//    cout << "dot(A,B)/(magA*magB)" << dot(A,B)/(magA*magB) << endl;
//    cout << "acos(dot(A,B)/(magA*magB)) " << SafeAcos(dot(A,B)/(magA*magB)) << endl;
    theta = M_PI - 0.5*SafeAcos(dot(A,B)/(magA*magB));
    phi = safeatan2(B(0),B(1));
    if (phi>=0 && phi<=M_PI){
      zeta = phi-M_PI;
    }else{
      zeta = phi+M_PI;
    }
    gamma = theta+zeta;
//    cout << "B(0)/magB" << B(0)/magB << endl;
//    cout << "theta" << theta << endl;
//    cout << "gamma" << gamma << endl;
    F(0) = sin(gamma);
    F(1) = cos(gamma);
    magF = dot(F,F);
    FdotR = dot(F,R);
//    cout << "F " << F(0) << " " << F(1) << endl;
//    cout << "FdotR " << FdotR << endl;
//    cout << "magF " << magF << endl;
//    cout << "F*elementii.adhesionmagnitude/magF " << F(0)*elementii.adhesionmagnitude/magF << F(1)*elementii.adhesionmagnitude/magF << endl;
    if (FdotR<0){
      elementii.fb += ((-1.0)*F*elementii.adhesionmagnitude/magF);
    }else{
      elementii.fb += (F*elementii.adhesionmagnitude/magF);
    }
//    elementii.fb.print();
  }
}
