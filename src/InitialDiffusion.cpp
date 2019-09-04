//
//  InitialDiffusion.cpp
//  ShapeIBM
//
//  Created by <author> on 23/08/2019.
//
//

#include "InitialDiffusion.hpp"
#include "InitialShape.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

void InitialDiffusion(mat& AccumGrid,const float& conc,const float& D,const float& tdif_max,const float& len,const float& h,const float& zeta,const int& Ng,const float& hg,const int& shapeflag){

  mat ConcentrationGrid = mat(Ng,Ng,fill::zeros);
  mat InitialGrid       = mat(Ng,Ng,fill::zeros);
  mat DeltaGrid         = mat(Ng,Ng,fill::zeros);
  float tdif=0;

  cout << "Initialising diffusion" << endl;

  //dtdif = pow(Tissue.hg,2)/(2*D);
  float dtdif = 1/(4*D);

  // Create initial diffusion concentration with InitialShape function.
  for (int i = 0; i < Ng; i++) {
    for (int j = 0; j < Ng; j++) {
      if (InitialShape(i,j,len,h,zeta,Ng,hg,shapeflag)){
        InitialGrid(i,j) = conc;
      }else{
        InitialGrid(i,j) = 0;
      }
    }
  }

  ConcentrationGrid=InitialGrid;
  while (tdif<tdif_max){
    for (int i=1;i<Ng-1;i++){
      for (int j=1;j<Ng-1;j++){
        float d2xPhi = ConcentrationGrid((i-1),j)-2*ConcentrationGrid(i,j)+ConcentrationGrid((i+1),j);
        float d2yPhi = ConcentrationGrid(i,(j-1))-2*ConcentrationGrid(i,j)+ConcentrationGrid(i,(j+1));
        DeltaGrid(i,j) = D*dtdif*(d2xPhi+d2yPhi);
      }
    }
    ConcentrationGrid = ConcentrationGrid+DeltaGrid;
    AccumGrid = AccumGrid+DeltaGrid;
    for (int i=0;i<Ng;i++){
      for (int j=0;j<Ng;j++){
        if (InitialGrid(i,j) > 0){
          AccumGrid(i,j)=0;
        }else{
          ConcentrationGrid(i,j)=0;
        }
      }
    }
    tdif = tdif+dtdif;
    system("clear");
    cout << "Diffusion progress: " << tdif << "/" << tdif_max << endl;
  }
}
