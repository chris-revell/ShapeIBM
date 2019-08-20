//
//  Initialise.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 22/04/2019.
//
//

#include "Initialise.hpp"
#include "element.hpp"
#include <armadillo>
#include <vector>
#include <math.h>
#include <InitialShape.hpp>
#include <ReadParameters.hpp>

using namespace std;
using namespace arma;

void Initialise(vector<ofstream>& files,vector<element>& Elements,int& Nb,int& Ng,float& rho,float& mu,float& re,float& tension,float& adhesion,float& dt,float& t_max,float& t_output,int& realtimeplot,float& xmin,float& xmax,float& hg,arma::mat& sb,arma::cube& xg,arma::mat& sg,arma::cube& fg,arma::cube& vg,arma::cube& ug,arma::mat& xbglobal,arma::mat& ubglobal,arma::mat& fbglobal){
  // System parameters
  int elementCount = 0;
  float dxjj_sq,D,tdif_max,dtdif,h,len,zeta;
  float tdif=0;
  float cen=0;
  float dimensions=1;

  mat PositionGrid,InitialGrid,Dgrid,AccumGrid;
  vec dxjj = vec(2,fill::zeros);
  vec dxn0 = vec(2,fill::zeros);
  vec dxn1 = vec(2,fill::zeros);

  ReadParameters(files[0],Ng,rho,mu,len,h,zeta,re,tension,adhesion,D,tdif_max,dt,t_max,t_output,realtimeplot);

  xg           = cube(Ng+1,Ng+1,2,fill::zeros);
  sg           = mat(Ng+1,Ng+1,fill::zeros);
  fg           = cube(Ng+1,Ng+1,2,fill::zeros);
  vg           = cube(Ng+1,Ng+1,2,fill::zeros);
  ug           = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal     = mat(2,0,fill::zeros);
  ubglobal     = mat(2,0,fill::zeros);
  fbglobal     = mat(2,0,fill::zeros);
  PositionGrid = mat(Ng,Ng,fill::zeros);
  InitialGrid = mat(Ng,Ng,fill::zeros);
  Dgrid        = mat(Ng,Ng,fill::zeros);
  AccumGrid    = mat(Ng,Ng,fill::zeros);
  xmin         =cen-static_cast<float>(dimensions);
  xmax         =cen+static_cast<float>(dimensions);
  sb(0,0)      = 0.0; // One sink point at the centre of each fluid grid edge
  sb(1,0)      = xmin;
  sb(0,1)      = 0.0;
  sb(1,1)      = xmax;
  sb(0,2)      = xmin;
  sb(1,2)      = 0.0;
  sb(0,3)      = xmax;
  sb(1,3)      = 0.0;
  hg           = float(xmax-xmin)/float(Ng);

  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg+hg/2;
      xg(ii,jj,1)=xmin+jj*hg+hg/2;
    }
  }
  // Write grid positions to file
  for (int ii=0;ii<Ng+1;ii++){
    files[5] << xg.slice(0).row(ii);
    files[6] << xg.slice(1).row(ii);
  }
  files[5].close();
  files[6].close();



  cout << "Initialising diffusion" << endl;


  //dtdif = pow(Tissue.hg,2)/(2*D);
  dtdif = 1/(4*D);

  // Create initial diffusion concentration with InitialShape function.
  for (int i = 0; i < Ng; i++) {
    for (int j = 0; j < Ng; j++) {
      if (InitialShape(i,j,len,h,zeta,Ng,hg)){
        InitialGrid(i,j) = 1;
      }else{
        InitialGrid(i,j) = 0;
      }
    }
  }

  PositionGrid=InitialGrid;
  while (tdif<tdif_max){
    for (int i=1;i<Ng-1;i++){
      for (int j=1;j<Ng-1;j++){
        float d2xPhi = PositionGrid((i-1),j)-2*PositionGrid(i,j)+PositionGrid((i+1),j);
        float d2yPhi = PositionGrid(i,(j-1))-2*PositionGrid(i,j)+PositionGrid(i,(j+1));
        Dgrid(i,j) = D*dtdif*(d2xPhi+d2yPhi);
      }
    }
    PositionGrid = PositionGrid+Dgrid;
    AccumGrid = AccumGrid+Dgrid;
    for (int i=0;i<Ng;i++){
      for (int j=0;j<Ng;j++){
        if (InitialGrid(i,j) > 0){
          AccumGrid(i,j)=0;
        }else{
          PositionGrid(i,j)=0;
        }
      }
    }
    tdif = tdif+dtdif;
    system("clear");
    cout << "Diffusion progress: " << tdif << "/" << tdif_max << endl;
  }


  // Create elements in regions of accumulated effector
  for (int ii=0; ii<Ng; ii++){
    for (int jj=0; jj<Ng; jj++){
      if (AccumGrid(ii,jj) > 0){
        Elements.push_back(element(elementCount,xg(ii,jj,0),xg(ii,jj,1),AccumGrid(ii,jj)));
        elementCount++;
      }
    }
  }

  // Set boundary arrays according to new element count
  Nb = elementCount;

  xbglobal.resize(2,elementCount);
  ubglobal.resize(2,elementCount);
  fbglobal.resize(2,elementCount);


  // Find neighbouring elements
  for (int ii=0; ii<elementCount; ii++){
    vec distances = vec(elementCount,fill::zeros);
    for (int jj=0; jj<elementCount; jj++){
      if (ii==jj){distances(jj) = 5000000;}
      else{
        dxjj = Elements[ii].pos - Elements[jj].pos;
        dxjj_sq = dot(dxjj,dxjj);
        distances(jj) = dxjj_sq;
      }
      uvec indices = sort_index(distances);
      Elements[ii].n0 = indices[0];
      Elements[ii].n1 = indices[1];
    }
  }  

  // Check neighbour results for errors
  for (int ii=0; ii<elementCount; ii++){
    element& elementii = Elements[ii];
    element& n0 = Elements[elementii.n0];
    element& n1 = Elements[elementii.n1];
    if (n0.n0 == ii || n0.n1 == ii){}
    else{cout << "neighbour error" << endl;
    cout << ii << " " << n0.n0 << " " << n0.n1 << " " << n1.n0 << " " << n1.n1 << endl;}
    if (n1.n0 == ii || n1.n1 == ii){}
    else{cout << "neighbour error" << endl;
    cout << ii << " " << n0.n0 << " " << n0.n1 << " " << n1.n0 << " " << n1.n1 << endl;}
  }

}
