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
    element(const float& v0, const float& v1, const int& cell, const int& Totalb,const float& initialx, const float& initialy, const int& n1, const int& n2);
    arma::vec pos = arma::vec(2,arma::fill::zeros);              // Element position
    arma::vec fb  = arma::vec(2,arma::fill::zeros);              // Forces on element
    arma::vec ub  = arma::vec(2,arma::fill::zeros);              // Previous velocity of element
    int label;                  // Global element label in xbglobal and fbglobal arrays
    int parent;                 // Label of cell to which element belongs
    float adhesionmagnitude;    // Adhesion magnitude of element
    std::vector<int> neighbours;// Indices of spatial neighbour elements in Cells[parent].Elements vector
    ~element();
protected:

};


#endif /* element_hpp */
