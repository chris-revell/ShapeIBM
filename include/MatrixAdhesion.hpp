//
//  MatrixAdhesion.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#ifndef MatrixAdhesion_hpp
#define MatrixAdhesion_hpp

#include <armadillo>
#include <vector>
#include <math.h>
#include "element.hpp"

void MatrixAdhesion(std::vector<element>& Elements,const int& Nb,const float& adhesion,const float& hg, const float& re);


#endif /* MatrixAdhesion_hpp */
