//
//  InitialDiffusion.hpp
//  ShapeIBM
//
//  Created by <author> on 23/08/2019.
//
//

#ifndef InitialDiffusion_hpp
#define InitialDiffusion_hpp

#include "InitialShape.hpp"
#include <armadillo>

void InitialDiffusion(arma::mat& AccumGrid,const float& conc,const float& D,const float& tdif_max,const float& len,const float& h,const float& zeta,const int& Ng,const float& hg,const int& shapeflag);

#endif /* InitialDiffusion_hpp */
