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
#include <tissue.hpp>
#include <math.h>
#include "ReadParameters.hpp"
#include "tissue.hpp"

using namespace std;
using namespace arma;

tissue Initialise(int& plotflag,float& t_m,float& t_out,ofstream& file){
  // System parameters
  int   Numg,Nb;
  int elementCount = 0;
  float dims,cen,Src,rho,mu,len,dt,tension,adhesion,dxjj_sq,h,A,alpha,D,tdif_max,dtdif,re;
  float tdif=0;
  mat PositionGrid,InitialGrid,Dgrid,AccumGrid;
  vec dxjj = vec(2,fill::zeros);
  vec dxn0 = vec(2,fill::zeros);
  vec dxn1 = vec(2,fill::zeros);

  ReadParameters(file,Numg,Nb,dims,cen,Src,rho,mu,len,dt,t_m,t_out,tension,adhesion,plotflag,h,A,alpha,D,tdif_max,re);

  // Create whole tissue system
  tissue Tissue = tissue(Numg,dims,Src,rho,mu,dt,cen,adhesion,tension,h,re);

  cout << "Initialising diffusion" << endl;

  float b = 2*A/((alpha+1)*h);
  //dtdif = pow(Tissue.hg,2)/(2*D);
  dtdif = 1/(4*D);
  PositionGrid = mat(Numg,Numg,fill::zeros);
  Dgrid        = mat(Numg,Numg,fill::zeros);
  AccumGrid    = mat(Numg,Numg,fill::zeros);

  for (int i = 0; i < Numg; i++) {
    float x = (i-Numg/2)*Tissue.hg;
    for (int j = 0; j < Numg; j++) {
      float y = (j-Numg/2)*Tissue.hg;
      //float y_1 = b*(alpha-1)*fabs(x)/(2*h)+b/2;
      //float y_0 = b*(1-alpha)*fabs(x)/(2*h)-b/2;
      //if (x>-h && x<h && y<y_1 && y>y_0){
      //  PositionGrid(i,j)=1;
      //}
      //else{
      //  PositionGrid(i,j)=0;
      //}
      float y_1 = b*(alpha-1)*x/(2*h)+b*(1+alpha)/4;
      float y_0 = b*(1-alpha)*x/(2*h)-b*(1+alpha)/4;
      if (x>(-h/2) && x<(h/2) && y<y_1 && y>y_0){
        PositionGrid(i,j)=1;
      }
      else{
        PositionGrid(i,j)=0;
      }
    }
  }
  InitialGrid = PositionGrid;
  while (tdif<tdif_max){
    for (int i=1;i<Numg-1;i++){
      for (int j=1;j<Numg-1;j++){
        float d2xPhi = PositionGrid((i-1),j)-2*PositionGrid(i,j)+PositionGrid((i+1),j);
        float d2yPhi = PositionGrid(i,(j-1))-2*PositionGrid(i,j)+PositionGrid(i,(j+1));
        Dgrid(i,j) = D*dtdif*(d2xPhi+d2yPhi);
      }
    }
    PositionGrid = PositionGrid+Dgrid;
    AccumGrid = AccumGrid+Dgrid;
    for (int i=0;i<Numg;i++){
      for (int j=0;j<Numg;j++){
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

  //PositionGrid.load("input/AccumGrid.txt",raw_ascii);
  //PGdim = PositionGrid.n_cols;

  for (int ii=0; ii<Numg; ii++){
    for (int jj=0; jj<Numg; jj++){
      if (AccumGrid(ii,jj) > 0){
        Tissue.Elements.push_back(element(elementCount,Tissue.xg(ii,jj,0),Tissue.xg(ii,jj,1),AccumGrid(ii,jj)));
        elementCount++;
      }
    }
  }

  Tissue.Nb = elementCount;
  Tissue.xbglobal.resize(2,elementCount);
  Tissue.ubglobal.resize(2,elementCount);
  Tissue.fbglobal.resize(2,elementCount);

  for (int ii=0; ii<elementCount; ii++){
    vec distances = vec(elementCount,fill::zeros);

    for (int jj=0; jj<elementCount; jj++){
      if (ii==jj){distances(jj) = 5000000;}
      else{
        dxjj = Tissue.Elements[ii].pos - Tissue.Elements[jj].pos;
        dxjj_sq = dot(dxjj,dxjj);
        distances(jj) = dxjj_sq;
      }
      uvec indices = sort_index(distances);
      Tissue.Elements[ii].n0 = indices[0];
      Tissue.Elements[ii].n1 = indices[1];
    }
  }

  for (int ii=0; ii<elementCount; ii++){
    element& elementii = Tissue.Elements[ii];
    element& n0 = Tissue.Elements[elementii.n0];
    element& n1 = Tissue.Elements[elementii.n1];
    if (n0.n0 == ii || n0.n1 == ii){}
    else{cout << "neighbour error" << endl;
    cout << ii << " " << n0.n0 << " " << n0.n1 << " " << n1.n0 << " " << n1.n1 << endl;}
    if (n1.n0 == ii || n1.n1 == ii){}
    else{cout << "neighbour error" << endl;
    cout << ii << " " << n0.n0 << " " << n0.n1 << " " << n1.n0 << " " << n1.n1 << endl;}
  }
  return Tissue;
}
