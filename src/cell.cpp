//
// Created by Christopher Revell on 24/10/2018.
//

#include "cell.hpp"
#include <armadillo>
#include <cmath>
#include "element.hpp"

using namespace arma;
using namespace std;

cell::cell(const int& cellnum, const int& Totalb, const int& NumBounds, const float& radius, const float& initialx, const float& initialy,const float& mesh,const float& tension) { // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
  Nb = NumBounds;
  label = cellnum;
  hb=2*M_PI/static_cast<float>(Nb); // Typical angular spacing between boundary elements given typical radius len
  com= vec(2,fill::zeros);    // Cell centre of mass
  ctension = tension; // Spring constant of boundary forces
  len=radius;
  hg=mesh;
  e=0.75;
  // Loop over initial element angles to set cartesian coordinates of all boundary points. Add constant value to set initial cell position.
  for (int ii=0;ii<Nb;ii++){
    r = len/sqrt(1-pow(e*cos(ii*hb),2));
    //r=len;
    Elements.push_back(element(0,0,label,Totalb+ii,r*cos(ii*hb)+initialx,r*sin(ii*hb)+initialy,((ii-1)%Nb+Nb)%Nb,((ii+1)%Nb+Nb)%Nb));
  }
  com(0) = initialx;
  com(1) = initialy;
}

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

float cell::CalculateVolume(){
  vec r1 = vec(3,fill::zeros);
  vec r2 = vec(3,fill::zeros);
  vec r3 = vec(3,fill::zeros);
  float volume = 0;
  for(int ii=0;ii<Nb;ii++){
    r1(0) = Elements[ii].pos(0)-com(0);
    r1(1) = Elements[ii].pos(1)-com(1);
    r2(0) = Elements[(ii+1)%Nb].pos(0)-com(0);
    r2(1) = Elements[(ii+1)%Nb].pos(1)-com(1);
    r3 = cross(r1,r2);
    volume = volume + 0.5*dot(r3,r3);
  }
  return volume;
}

cell::~cell() {}
