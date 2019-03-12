//
//  AdjacentForces.hpp
//  ImmersedBoundary
//
//  Created by <author> on 20/02/2019.
//
//

#ifndef AdjacentForces_hpp
#define AdjacentForces_hpp

#include <stdio.h>
#include <armadillo>
#include <cmath>
#include "element.hpp"
#include "cell.hpp"

void AdjacentForces(cell& Cell, const float& time, const float& diffusionconstant);


#endif /* AdjacentForces_hpp */
