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
    //arma::vec ub  = arma::vec(2,arma::fill::zeros);              // Velocity of element
    //int label;                  // Global element label in xbglobal and fbglobal arrays
    //int parent;                 // Label of cell to which element belongs
    //float adhesionmagnitude;    // Adhesion magnitude of element
    //float baselineadhesion;     // Baseline adhesion for cell
    //float normalisationfactor;  // Factor to normalise adhesion by area
    //std::vector<int> neighbours;// Indices of spatial neighbour elements in Cells[parent].Elements vector
    ~element();
protected:

};


#endif /* element_hpp */
