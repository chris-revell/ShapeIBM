//
//  ReadParameters.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 20/03/2019.
//
//

#ifndef ReadParameters_hpp
#define ReadParameters_hpp

#include <stdio.h>
#include "ReadParameters.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void ReadParameters(int& Numg,int& Nb,float& dims,float& cen,float& Src,float& rho,float& mu,float& len,int& Numcells,float& dt,float& t_max,float& t_output,float& tension,float& adhesion,int& realtimeplot);

#endif /* ReadParameters_hpp */
