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
#include <armadillo>
#include "element.hpp"
#include <vector>
#include <math.h>

void MatrixAdhesion(std::vector<element>& Elements,const int& Nb,const float& adhesion,const float& hg);


#endif /* MatrixAdhesion_hpp */
