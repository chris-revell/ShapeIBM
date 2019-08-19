//
//  Initialise.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 22/04/2019.
//
//

#ifndef Initialise_hpp
#define Initialise_hpp


#include "element.hpp"
#include <armadillo>
#include <vector>
#include <tissue.hpp>
#include <math.h>

tissue Initialise(int& plotflag,float& t_m,float& t_out,std::ofstream& file);

#endif /* Initialise_hpp */
