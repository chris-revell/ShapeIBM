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
#include "element.hpp"
#include "smallfunctions.hpp"

class tissue {
private:

public:
    tissue(const int& GridSize,const int& dimensions,const int& boundarypoints,const float& sourcestrength, const float& density, const float& viscocity,const float& timestep); // Constructor
    void AddCell(const float& len, const float& initialx, const float& initialy,float& tension); // Function to add a cell object to the tissue
    void CombineBoundaries(void);
    void UpdateSources(void);
    void BoundaryRefinement(void);
    void UpdatePositions(void);
    void MatrixAdhesions(void);
    arma::cube xg;                // Fluid grid
    arma::cube xNb;               // Number of boundary elements within each fluid grid component
    arma::mat  sg;                // Grid source distribution
    arma::cube fg;                // Fluid forces
    arma::cube vg;                // Fluid velocities
    arma::cube ug;                // Previous fluid velocities
    arma::mat  sb;                // Sink and source positions
    arma::mat  sbb;               // Sink and source magnitudes
    arma::mat  xbglobal;          // Positions of all boundary points in the system
    arma::mat  fbglobal;          // Forces on all boundary points in the system
    arma::mat  ubglobal;          // Velocities of all boundary points in the system
    arma::mat  indices;           // Velocities of all boundary points in the system
    int Ng;                       // Fluid grid dimensions of the system
    int Nc;                       // Number of cells in the system
    int Nbcell;                   // Number of boundary points per cell
    int Nb;                       // Total number of boundary points in system
    int Nbs;                      // Number of sinks and sources in the system = Nc+4
    std::vector<cell> Cells;      // Vector of cell objects in system
    float xmin;                   // Fluid domain dimensions
    float xmax;                   // Fluid domain dimensions
    float hg;                     // Fluid mesh width
    float Src;                    // Source strength/cell growth rate
    float dt;                     // Time interval between steps
    float rho;                    // Fluid density
    float mu;                     // Fluid drag factor
    ~tissue();                    // Destructor
protected:

};


#endif /* tissue_hpp */
