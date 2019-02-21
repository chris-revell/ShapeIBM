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
void AdjacentForces(cell& Cell) {
  float dx1,dy1,dx2,dy2,r1,r2;
  for (int ii=0;ii<Cell.Nb;ii++){
    // Access the labels of the neighbours for element ii and extract the corresponding element positions.
    // Use these to find the separation of ii from its neighbours in the x and y directions.
    dx1 = Cell.Elements[Cell.Elements[ii].neighbours[0]].pos(0)-Cell.Elements[ii].pos(0);
    dy1 = Cell.Elements[Cell.Elements[ii].neighbours[0]].pos(1)-Cell.Elements[ii].pos(1);
    dx2 = Cell.Elements[Cell.Elements[ii].neighbours[1]].pos(0)-Cell.Elements[ii].pos(0);
    dy2 = Cell.Elements[Cell.Elements[ii].neighbours[1]].pos(1)-Cell.Elements[ii].pos(1);
    // Find separation distances from x and y values.
    r1=sqrt(pow(dx1,2)+pow(dy1,2));
    r2=sqrt(pow(dx2,2)+pow(dy2,2));
    // Calculate forces on element ii from
    //Elements[ii].fb(0) = ctension*((r1-hb)*dx1/r1+(r2-hb)*dx2/r2);
    //Elements[ii].fb(1) = ctension*((r1-hb)*dy1/r1+(r2-hb)*dy2/r2);
    Cell.Elements[ii].fb(0) = Cell.ctension*(dx1/r1+dx2/r2);
    Cell.Elements[ii].fb(1) = Cell.ctension*(dy1/r1+dy2/r2);
  }
}
