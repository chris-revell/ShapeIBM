//
//  BoundToGrid2.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//
#ifndef BOUNDTOGRID2_H
#define BOUNDTOGRID2_H

#include <armadillo>
#include "smallfunctions.hpp"

void BoundToGrid2(arma::cube& fg,const int& Nb,const arma::mat& fbglobal,const arma::mat& xbglobal,const float& xmin,const float& xmax,const float& hg,const int& Ng);


#endif
