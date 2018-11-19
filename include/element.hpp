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
    arma::vec pos;              // Element position
    arma::vec internalforce;    // Intra-cell forces on element
    arma::vec ub;               // Previous velocity of element
    int label;                  // Global element label in xbglobal and fbglobal arrays
    int parent;                 // Label of cell to which element belongs
    std::vector<int> neighbours;// Indices of spatial neighbour elements in Cells[parent].Elements vector
    ~element();
protected:

};


#endif /* element_hpp */
