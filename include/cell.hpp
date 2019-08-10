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
  cell(const int& cellnum, const int& Totalb, const int& NumBounds, const float& radius, const float& initialx, const float& initialy,const float& mesh,const float& tension,const float& adhesion); // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
  void AdjacentForces(void);
  void OppositeForces(void);
  void UpdateCom(void);
  void NormaliseAdhesion(void);
  std::vector<element> Elements; // Vector containing all element objects within this cell object. *Includes deleted elements*
  //std::vector<int> ElementLabels;// Vector containing all element objects *that currently exist*
  arma::vec com;                 // Cell centre of mass
  float ctension;                // Spring constant of boundary forces
  float adhesionmagnitude;       // Baseline adhesion magnitude
  float hb;                      // Typical angular spacing between boundary elements given typical radius len
  float hg;                      // Global fluid mesh spacing
  float len;                     // Typical cell radius
  int Nb;                        // Number of boundary points in cell
  int NbT;                       // Total length of the Elements array including those that have been deleted.
  int label;                     // Identifying label of cell
  float CalculateVolume();
  ~cell();                       //Destructor
protected:

};


#endif /* test_hpp */
