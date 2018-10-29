//
// Created by Christopher Revell on 24/10/2018.
// Class to contain data for a cell object in the IBM implementation, including boundary element positions,

#ifndef CELL_H
#define CELL_H

#include <armadillo>

class cell {
private:
public:
    cell(const int& NumBounds, const float& len, const float& initialx,const float& initialy); // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
    arma::mat xb; // Positions of all boundary points in cell
    arma::mat fb; // Forces on all boundary points in cell arising from interactions with other boundary points
    float hb;     // Typical angular spacing between boundary elements given typical radius len
    int Nb;       // Number of boundary points in cell
    void AdjacentForces(void);
    float corticaltension; // Spring constant of boundary forces
    ~cell(); //Destructor
protected:

};


#endif /* test_hpp */
