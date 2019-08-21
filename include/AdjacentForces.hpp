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
#include <math.h>
#include "element.hpp"
#include <vector>

void AdjacentForces(std::vector<element>& Elements,const float& re,const int& Nb,const float& tension);


#endif /* AdjacentForces_hpp */
