//
//  tissue.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 26/10/2018.
//
//

#ifndef tissue_hpp
#define tissue_hpp
#include <vector>
#include <cell>

class tissue {
private:
public:
    tissue(const int& Ng);
    arma::mat xbglobal; // Positions of all boundary points in the system
    arma::mat sources;  // Positions of all fluid sources in the system
    arma::mat sinks;    // Positions of all fluid sinks in the system
    int Nc;             // Number of cells in the system
    vector<cell> Cells; // Vector of cell objects in system
    ~tissue();
protected:

};


#endif /* tissue_hpp */
