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
    tissue(const int& GridSize,const int& dimensions,const int& boundarypoints,const float& sourcestrength); // Constructor
    void AddCell(const float& len, const float& initialx, const float& initialy); // Function to add a cell object to the tissue
    void CombineBoundaries(void);
    void UpdateSources(void);
    arma::cube xg;           // Fluid grid
    arma::mat  sg;
    arma::cube fg;           // Fluid forces
    arma::cube vg;           // Fluid velocities
    arma::cube ug;           // Previous fluid velocities
    arma::mat  sb;           // Sink and source positions
    arma::mat  sbb;          // Sink and source magnitudes
    arma::mat  xbglobal;     // Positions of all boundary points in the system
    arma::mat  fbglobal;     // Forces on all boundary points in the system
    arma::mat  ubglobal;     // Velocities of all boundary points in the system
    int Ng;                  // Fluid grid dimensions of the system
    int Nc;                  // Number of cells in the system
    int Nbcell;              // Number of boundary points per cell
    int Nb;                  // Total number of boundary points in system
    int Nbs;                 // Number of sinks and sources in the system = Nc+4
    std::vector<cell> Cells; // Vector of cell objects in system
    float xmin;              // Fluid domain dimensions
    float xmax;              // Fluid domain dimensions
    float hg;                // Fluid mesh width
    float Src;               // Source strength/cell growth rate
    ~tissue();               // Destructor
protected:

};


#endif /* tissue_hpp */
