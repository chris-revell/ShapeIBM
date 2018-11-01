//
//  GridToBound.hpp
//  ImmersedBoundary
//
//  Created by <author> on 26/10/2018.
//
//
#ifndef GRIDTOBOUND_H
#define GRIDTOBOUND_H

#include <armadillo>
void GridToBound(arma::mat& fb,const arma::mat& xb,const int& Nb,const arma::cube& fg,const int& Ng,const float& hdl,const float& hg,const float& xmn,const float& xmx);


#endif
