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
    element(const float& initialx, const float& initialy, const int& n1, const int& n2);
    arma::vec pos;  // Element position
    int label;      // Unique element label
    int parent;     // Label of cell to which element belongs
    std::vector<int> neighbours; // Labels of elements with with this element shares a tension interaction in the boundary
    ~element();
protected:

};


#endif /* element_hpp */
