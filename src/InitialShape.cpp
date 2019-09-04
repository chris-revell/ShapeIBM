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

bool InitialShape(const int& i,const int& j,const float& len, const float& h, const float& zeta,const int& Ng,const float& hg,const int& shapeflag){

  if (shapeflag==0){
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
  }

  else{
    //float A = h*len;
    //float b = 2*A/((zeta+1)*len);
    float x = (i-Ng/2)*hg;
    float y = (j-Ng/2)*hg;
    float y_1 = h*(zeta-1)*x/(2*len)+h*(1+zeta)/4;
    float y_0 = h*(1-zeta)*x/(2*len)-h*(1+zeta)/4;
    if (x>(-len) && x<(len) && y<y_1 && y>y_0){
      return true;
    }
    else{
      return false;
    }
  }








}
