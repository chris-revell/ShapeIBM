//
//  LocalToGlobal.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#ifndef LocalToGlobal_hpp
#define LocalToGlobal_hpp

#include "element.hpp"
#include <armadillo>
#include <vector>

void LocalToGlobal(const std::vector<element>& Elements,arma::mat& xbglobal,arma::mat& fbglobal,arma::mat& ubglobal,const int& Nb);

#endif /* LocalToGlobal_hpp */
