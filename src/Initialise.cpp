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
#include <ReadParameters.hpp>
#include <InitialDiffusion.hpp>

using namespace std;
using namespace arma;

void Initialise(vector<ofstream>& files,vector<element>& Elements,int& Nb,int& Ng,float& rho,float& mu,float& re,float& tension,float& adhesion,float& dt,float& t_max,float& t_output,int& realtimeplot,float& xmin,float& xmax,float& hg,arma::cube& xg,arma::mat& sg,arma::cube& fg,arma::cube& vg,arma::cube& ug,arma::mat& xbglobal,arma::mat& ubglobal,arma::mat& fbglobal){
  // System parameters
  int elementCount = 0;
  float dxjj_sq,D,tdif_max,h,len,zeta,conc;
  float cen=0;
  float dimensions=1;

  mat AccumGrid;
  vec dxjj = vec(2,fill::zeros);
  vec dxn0 = vec(2,fill::zeros);
  vec dxn1 = vec(2,fill::zeros);

  ReadParameters(files[0],Ng,rho,mu,len,h,zeta,re,tension,adhesion,D,conc,tdif_max,dt,t_max,t_output,realtimeplot);

  xg           = cube(Ng+1,Ng+1,2,fill::zeros);
  sg           = mat(Ng+1,Ng+1,fill::zeros);
  fg           = cube(Ng+1,Ng+1,2,fill::zeros);
  vg           = cube(Ng+1,Ng+1,2,fill::zeros);
  ug           = cube(Ng+1,Ng+1,2,fill::zeros);
  xbglobal     = mat(2,0,fill::zeros);
  ubglobal     = mat(2,0,fill::zeros);
  fbglobal     = mat(2,0,fill::zeros);
  AccumGrid    = mat(Ng,Ng,fill::zeros);
  xmin         =cen-static_cast<float>(dimensions);
  xmax         =cen+static_cast<float>(dimensions);
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

  InitialDiffusion(AccumGrid,conc,D,tdif_max,len,h,zeta,Ng,hg);

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
