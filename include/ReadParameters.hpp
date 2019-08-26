//
//  ReadParameters.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 20/03/2019.
//
//

#ifndef ReadParameters_hpp
#define ReadParameters_hpp

#include "ReadParameters.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void ReadParameters(std::ofstream& file,int& Ng,float& rho,float& mu,float& len,float& h,float& zeta,float& re,float& tension,float& adhesion,float& D,float& conc,float& tdif_max,float& dt,float& t_max,float& t_output,int& realtimeplot);

#endif /* ReadParameters_hpp */
