//
// Created by Christopher Revell on 24/10/2018.
//

#include "cell.hpp"
#include <armadillo>
#include <cmath>
#include "element.hpp"

using namespace arma;
using namespace std;

cell::cell(const int& cellnum, const int& Totalb, const int& NumBounds, const float& radius, const float& initialx, const float& initialy,const float& mesh,const float& tension,const float& adhesion) { // Constructor takes number of boundary points, typical radius, and x,y positions of centre of mass
  Nb = NumBounds;
  NbT = Nb;
  label = cellnum;
  hb=2*M_PI/static_cast<float>(Nb); // Typical angular spacing between boundary elements given typical radius len
  com= vec(2,fill::zeros);          // Cell centre of mass
  ctension = tension;               // Spring constant of boundary forces
  len=radius;
  hg=mesh;
  e=0.9;
  adhesionmagnitude = adhesion;
  // Loop over initial element angles to set cartesian coordinates of all boundary points. Add constant value to set initial cell position.
  for (int ii=0;ii<Nb;ii++){
    r = len/sqrt(1-pow(e*cos(ii*hb),2));
    Elements.push_back(element(label,ii,r*cos(ii*hb)+initialx,r*sin(ii*hb)+initialy,((ii-1)%Nb+Nb)%Nb,((ii+1)%Nb+Nb)%Nb,adhesion));
  }
  com(0) = initialx;
  com(1) = initialy;
}

// Function to update cell centre of mass from boundary positions
void cell::UpdateCom(){
  vec sum = vec(2,fill::zeros);
  for(int ii=0;ii<Nb;ii++){
    sum = sum + Elements[ii].pos;
  }
  com = sum/Nb;
}

float cell::CalculateVolume(){
  vec r1 = vec(3,fill::zeros);
  vec r2 = vec(3,fill::zeros);
  vec r3 = vec(3,fill::zeros);
  float volume = 0;
  for(int ii=0;ii<Nb;ii++){
    element& elementii = Elements[ii];
    element& n1 = Elements[elementii.neighbours[1]];
    r1(0) = elementii.pos(0)-com(0);
    r1(1) = elementii.pos(1)-com(1);
    //r1 = Elements[ElementLabels[ii]].pos-com;
    r2(0) = n1.pos(0)-com(0);
    r2(1) = n1.pos(1)-com(1);
    //r2 = Elements[ElementLabels[(ii+1)%Nb]].pos-com;
    r3 = cross(r1,r2);
    volume = volume + 0.5*sqrt(dot(r3,r3));
  }
  return volume;
}

void cell::NormaliseAdhesion(void){
  for (int ii=0;ii<Nb;ii++){
    element& elementii = Elements[ii];
    element& n0 = Elements[elementii.neighbours[0]];
    element& n1 = Elements[elementii.neighbours[1]];
    vec d0 = elementii.pos-n0.pos;
    vec d1 = elementii.pos-n1.pos;
    float r0 = sqrt(dot(d0,d0));
    float r1 = sqrt(dot(d1,d1));
    elementii.normalisationfactor = (r0+r1)/2.0;
  }
}

cell::~cell() {}
