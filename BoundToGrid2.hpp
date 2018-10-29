//
//  BoundToGrid2.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#include <armadillo>
void BoundToGrid2(arma::cube& sg,const arma::mat& xb,const arma::mat& sb,const int& Nb,const int& Ng,const float& hdl,const float& hg,const float& hb,const float& xmin,const float& xmax);
