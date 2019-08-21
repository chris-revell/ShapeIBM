//
//  GlobalToLocal.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 01/11/2018.
//
//
#ifndef GLOBALTOLOCAL_H
#define GLOBALTOLOCAL_H

#include "GlobalToLocal.hpp"
#include <armadillo>
#include <vector>

void GlobalToLocal(std::vector<element>& Elements,const arma::mat& xbglobal,const int& Nb);

#endif
