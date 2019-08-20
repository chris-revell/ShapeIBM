//
//  element.hpp
//  ImmersedBoundary
//
//  Created by <author> on 07/11/2018.
//
//

#ifndef element_hpp
#define element_hpp

#include <armadillo>
#include <vector>

class element {
private:
public:
    element(const int& Label, const float& initialx, const float& initialy, const float& accum);
    void SetAdhesion(void);
    arma::vec pos        = arma::vec(2,arma::fill::zeros);       // Element position
    arma::vec initialpos = arma::vec(2,arma::fill::zeros);       // Initial element position
    arma::vec fb         = arma::vec(2,arma::fill::zeros);       // Forces on element
    int n0;
    int n1;
    float accumulatedEffector;    
    ~element();
protected:

};


#endif /* element_hpp */
