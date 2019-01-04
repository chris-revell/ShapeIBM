//
// Created by Christopher Revell on 24/10/2018.
// Class to contain data for a cell object in the IBM implementation, including boundary element positions,

#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include "element.hpp"

class cell {
private:
  float r;
  float e;
public:
  cell(const int& cellnum, const int& Totalb, const int& NumBounds, const float& radius, const float& initialx, const float& initialy,const float& mesh,const float& tension); // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
  void AdjacentForces(void);
  void OppositeForces(void);
  void UpdateCom(void);
  std::vector<element> Elements; // Vector containing all element objects within this cell object
  arma::mat xb;                  // Positions of all boundary points in cell
  arma::mat fb;                  // Forces on all boundary points in cell arising from interactions with other boundary points
  arma::vec com;                 // Cell centre of mass
  float ctension;         // Spring constant of boundary forces
  float hb;                      // Typical angular spacing between boundary elements given typical radius len
  float hg;                      // Global fluid mesh spacing
  float len;                     // Typical cell radius
  int Nb;                        // Number of boundary points in cell
  int label;                     // Identifying label of cell
  float CalculateVolume();
  ~cell();                       //Destructor
protected:

};


#endif /* test_hpp */
