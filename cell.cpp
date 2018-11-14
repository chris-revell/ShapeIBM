//
// Created by Christopher Revell on 24/10/2018.
//

#include "cell.hpp"
#include <armadillo>
#include <cmath>
#include "element.hpp"

using namespace arma;
using namespace std;

cell::cell(const int& cellnum, const int& Totalb, const int& NumBounds, const float& radius, const float& initialx, const float& initialy,const float& mesh) { // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
  Nb = NumBounds;
  label = cellnum;
  hb=2*M_PI/static_cast<float>(Nb); // Typical angular spacing between boundary elements given typical radius len
  xb = mat(2,Nb,arma::fill::zeros); // Positions of all boundary points in cell
  fb = mat(2,Nb,arma::fill::zeros); // Forces on all boundary points in cell arising from interactions with other boundary points
  com= vec(2,arma::fill::zeros);    // Cell centre of mass
  corticaltension = 0.01; // Spring constant of boundary forces
  len=radius;
  hg=mesh;
  // Loop over initial element angles to set cartesian coordinates of all boundary points. Add constant value to set initial cell position.
  for (int ii=0;ii<Nb;ii++){
    Elements.push_back(element(label,Totalb+ii,len*cos(ii*hb)+initialx,len*sin(ii*hb)+initialy,((ii-1)%Nb+Nb)%Nb,((ii+1)%Nb+Nb)%Nb));
  }
  com(0) = initialx;
  com(1) = initialy;
}

// Function to calculate forces between neighbouring boundary elements within a cell.
void cell::AdjacentForces() {
  float dx1,dy1,dx2,dy2,r1,r2;
  if (xb.n_cols < Nb){
    xb.resize(2,Nb);
    fb.resize(2,Nb);
  }
  for (int ii=0;ii<Nb;ii++){
    // Access the labels of the neighbours for element ii and extract the corresponding element positions.
    // Use these to find the separation of ii from its neighbours in the x and y directions.
    dx1 = Elements[Elements[ii].neighbours[0]].pos(0)-Elements[ii].pos(0);
    dy1 = Elements[Elements[ii].neighbours[0]].pos(1)-Elements[ii].pos(1);
    dx2 = Elements[Elements[ii].neighbours[1]].pos(0)-Elements[ii].pos(0);
    dy2 = Elements[Elements[ii].neighbours[1]].pos(1)-Elements[ii].pos(1);
    // Find separation distances from x and y values.
    r1=sqrt(pow(dx1,2)+pow(dy1,2));
    r2=sqrt(pow(dx2,2)+pow(dy2,2));
    // Calculate forces on element ii from
    fb(0,ii)=corticaltension*((r1-hb)*dx1/r1+(r2-hb)*dx2/r2);
    fb(1,ii)=corticaltension*((r1-hb)*dy1/r1+(r2-hb)*dy2/r2);
  }
}

//--------------------------------------------------------------------//
// define opposite forces between opposite boundary points            //
// Again Hookean springs with spring constant corticaltension and this time       //
// equilibrium radius Lrest=2*len where len is the tyM_PIcal cell radius//
//--------------------------------------------------------------------//
//void cell::OppositeForces(){
//  float Lrest;
//  float dx1;
//  float dy1;
//  float r1;
//  int Nb2;
//
//  Lrest=2*len;
//  Nb2=floor(Nb/4);
//  for (int ii=0; ii<Nb2; ii++){
//    dx1=xb(0,3*Nb2-ii)-xb(0,ii);
//    dy1=xb(1,3*Nb2-ii)-xb(1,ii);
//    r1 = sqrt(pow(dx1,2)+pow(dy1,2));
//    fb(0,ii)       = fb(0,ii)      +corticaltension*abs(r1-Lrest)*dx1/r1;
//    fb(1,ii)       = fb(1,ii)      +corticaltension*abs(r1-Lrest)*dy1/r1;
//    fb(0,3*Nb2-ii) = fb(0,3*Nb2-ii)-corticaltension*abs(r1-Lrest)*dx1/r1;
//    fb(1,3*Nb2-ii) = fb(1,3*Nb2-ii)-corticaltension*abs(r1-Lrest)*dy1/r1;
//  }
//}// function OppositeForces
//--------------------------------------------------------------------//

// Function to update cell centre of mass from boundary positions
void cell::UpdateCom(){
  float xsum=0;
  float ysum=0;
  for(int ii=0;ii<Nb;ii++){
    xsum = xsum + Elements[ii].pos(0);
    ysum = ysum + Elements[ii].pos(1);
  }
  com(0) = xsum/Nb;
  com(1) = ysum/Nb;
}





cell::~cell() {}
