//
// Created by Christopher Revell on 24/10/2018.
// Class to contain data for a cell object in the IBM implementation, including boundary element positions,

#ifndef CELL_H
#define CELL_H

#include <armadillo>

class cell {
private:
public:
    cell(const int& Nb, const float& len, const float& initialpositionx,const float& initialpositiony);
    arma::mat xb;
    ~cell();
protected:

};


#endif /* test_hpp */
