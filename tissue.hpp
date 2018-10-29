//
//  tissue.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#ifndef tissue_hpp
#define tissue_hpp

#include "tissue.hpp"
#include <armadillo>
#include <vector>
#include "cell.hpp"
#include "smallfunctions.hpp"

class tissue {
private:

public:
    tissue(const int& GridSize,const int& dimensions);
    void AddCell(const int& BoundPoints, const float& len, const float& initialx, const float& initialy);
    arma::cube xg;           // Fluid grid
    arma::mat  sg;           // Fluid grid
    arma::mat xbglobal;      // Positions of all boundary points in the system
    //arma::mat sources;     // Positions of all fluid sources in the system
    //arma::mat sinks;       // Positions of all fluid sinks in the system
    int Ng;                  // Fluid grid dimensions of the system
    int Nc;                  // Number of cells in the system
    int Nb;                  // Number of boundary points per cell
    std::vector<cell> Cells; // Vector of cell objects in system
    float xmin;              // Fluid domain dimensions
    float xmax;              // Fluid domain dimensions
    float hg;                // Fluid mesh width
    ~tissue();
protected:

};


#endif /* tissue_hpp */
