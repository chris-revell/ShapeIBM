//
//  InitialShape.cpp
//  ShapeIBM
//
//  Created by <author> on 19/08/2019.
//
//

#include "InitialShape.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

bool InitialShape(const int& i,const int& j,const float& len, const float& h, const float& zeta,const int& Ng,const float& hg){
  float x,y,ylim;
  x = (i-Ng/2)*hg;
  y = (j-Ng/2)*hg;
  if (fabs(x)>len){
    return false;
  }else{
    if (x < -zeta*len){
      ylim = h*(len+x)/(len*(1-zeta));
    }else{
      ylim = h*(len-x)/(len*(1+zeta));
    }
    if (fabs(y)<ylim){
      return true;
    }else{
      return false;
    }
  }





  

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









}
