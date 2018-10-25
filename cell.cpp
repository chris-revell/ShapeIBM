//
// Created by Christopher Revell on 24/10/2018.
//

#include "cell.hpp"
#include "armadillo"
#include "cmath"

cell::cell(const int& Nb, const float& len, const float& initialpositionx, const float& initialpositiony) {

    float hb=2*M_PI/static_cast<float>(Nb);

    xb = arma::mat(2,Nb,arma::fill::zeros);

    for (int ii=0;ii<Nb;ii++){
        xb(0,ii) = len*cos(ii*hb)+initialpositionx;
        xb(1,ii) = len*sin(ii*hb)+initialpositiony;
    }
}
cell::~cell() {}
