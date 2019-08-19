//
//  AdjacentForces.cpp
//  ImmersedBoundary
//
//  Created by <author> on 20/02/2019.
//
//

#include "AdjacentForces.hpp"
#include <armadillo>
#include <math.h>
#include "tissue.hpp"
#include "element.hpp"
#include <vector>

using namespace std;
using namespace arma;

// Function to calculate forces between neighbouring boundary elements within a cell.
void AdjacentForces(tissue& Tissue) {
  float r0,r1,localSpringConstant0,localSpringConstant1;
  vec d0 = vec(2,fill::zeros);
  vec d1 = vec(2,fill::zeros);
  float& re = Tissue.r_e;

  for (int ii=0;ii<Tissue.Nb;ii++){
    element& elementii = Tissue.Elements[ii];
    element& n0 = Tissue.Elements[elementii.n0];
    element& n1 = Tissue.Elements[elementii.n1];
    // Access the labels of the neighbours for element ii and extract the corresponding element positions.
    // Use these to find the separation of ii from its neighbours in the x and y directions.
    d0 = n0.pos-elementii.pos;
    d1 = n1.pos-elementii.pos;
    // Find separation distances from x and y values.
    r0=sqrt(dot(d0,d0));
    r1=sqrt(dot(d1,d1));
    localSpringConstant0 = Tissue.ctension*(n0.accumulatedEffector+elementii.accumulatedEffector)/2.0;
    localSpringConstant1 = Tissue.ctension*(n1.accumulatedEffector+elementii.accumulatedEffector)/2.0;
    elementii.fb = localSpringConstant0*(r0-re)*d0/r0+localSpringConstant1*(r1-re)*d1/r1;
  }
}
