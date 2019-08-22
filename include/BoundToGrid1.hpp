//
//  BoundToGrid1.hpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//
#ifndef BOUNDTOGRID1_H
#define BOUNDTOGRID1_H

#include <armadillo>
#include "smallfunctions.hpp"

void BoundToGrid1(arma::mat& sg,const int& Nbs,const arma::mat& sb,const arma::mat& sbb,const float& xmin,const float& xmax,const float& hg,const int& Ng);


#endif
