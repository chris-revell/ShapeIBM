//
//  AdjacentForces.cpp
//  ImmersedBoundary
//
//  Created by <author> on 20/02/2019.
//
//

#include "AdjacentForces.hpp"
#include <armadillo>
#include <cmath>
#include "element.hpp"
#include "cell.hpp"

using namespace std;
using namespace arma;

// Function to calculate forces between neighbouring boundary elements within a cell.
void AdjacentForces(cell& Cell,const float& time) {
  float r0,r1,localtension,s0,s1,diffusiondistance;
  vec d0 = vec(2,fill::zeros);
  vec d1 = vec(2,fill::zeros);
  vector<element>& Elements  = Cell.Elements;
  vector<int>& ElementLabels = Cell.ElementLabels;
  for (int ii=0;ii<Cell.Nb;ii++){
    element& elementii = Elements[ElementLabels[ii]];
    element& n0 = Elements[elementii.neighbours[0]];
    element& n1 = Elements[elementii.neighbours[1]];
    // Access the labels of the neighbours for element ii and extract the corresponding element positions.
    // Use these to find the separation of ii from its neighbours in the x and y directions.
    d0 = n0.pos-elementii.pos;
    d1 = n1.pos-elementii.pos;
    // Find separation distances from x and y values.
    r0=sqrt(dot(d0,d0));
    r1=sqrt(dot(d1,d1));
    elementii.fb = Cell.ctension*(d0/r0+d1/r1);
  }
}
