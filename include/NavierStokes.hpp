//
//  NavierStokes.hpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//
#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H

#include <armadillo>
#include "smallfunctions.hpp"

void NavierStokes(arma::cube& ug,arma::cube& vg,const arma::cube& fg,const arma::mat& sg,const float& Ng,const float& mu,const float& rho,const float& hg,const float& dt);

#endif
