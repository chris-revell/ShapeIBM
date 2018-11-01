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

void NavierStokes(arma::cube& vg,const arma::cube& ug,const arma::cube& fg,const arma::mat& sg,const int& Ng,const float& rho,const float& mu,const float& dt,const float& hg);


#endif
