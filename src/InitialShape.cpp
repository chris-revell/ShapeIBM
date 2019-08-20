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
}
