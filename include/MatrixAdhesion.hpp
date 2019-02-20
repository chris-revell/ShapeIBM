//
//  MatrixAdhesion.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#ifndef MatrixAdhesion_hpp
#define MatrixAdhesion_hpp

#include <stdio.h>
#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "element.hpp"
#include <cmath>

double SafeAcos (double x);
double safeatan2 (double x,double y);
void MatrixAdhesion(cell& Cell);


#endif /* MatrixAdhesion_hpp */
