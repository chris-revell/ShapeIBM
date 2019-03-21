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

void MatrixAdhesion(tissue& Tissue, const float& diffusionconstant, const float& time);


#endif /* MatrixAdhesion_hpp */
