//
//  GridToBound.hpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//
#ifndef GRIDTOBOUND_H
#define GRIDTOBOUND_H

#include "tissue.hpp"
#include <armadillo>
#include "smallfunctions.hpp"

void GridToBound(arma::mat& ubglobal,const arma::mat& xbglobal,const arma::cube& vg,const float& hg,const float& xmin,const float& xmax,const int& Nb,const int& Ng);

#endif
