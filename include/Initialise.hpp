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
#include <math.h>
#include <ReadParameters.hpp>
#include <InitialDiffusion.hpp>

void Initialise(std::vector<std::ofstream>& files,std::vector<element>& Elements,int& Nb,int& Ng,float& rho,float& mu,float& re,float& tension,float& adhesion,float& dt,float& t_max,float& t_output,int& realtimeplot,float& xmin,float& xmax,float& hg,arma::cube& xg,arma::mat& sg,arma::cube& fg,arma::cube& vg,arma::cube& ug,arma::mat& xbglobal,arma::mat& ubglobal,arma::mat& fbglobal);

#endif /* Initialise_hpp */
