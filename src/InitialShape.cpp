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

  float x = (i-Ng/2)*hg;
  float y = (j-Ng/2)*hg;
  float y_1,y_0;

  if (shapeflag==0){

    if (x < zeta*2.0*len-len){
      y_1 = x*h/(zeta*2.0*len) + h/(2.0*zeta);
      y_0 = -y_1;
      if (y<y_1 && y>y_0){
        return true;
      }
      else{
        return false;
      }
    }else{
      y_1 = -x*h/(2.0*len*(1.0-zeta)) + h/(2.0*(1.0-zeta)) ;//(h/(zeta-1.0))*(x/len-2.0*h);
      y_0 = -y_1;//-(h/(zeta-1.0))*(x/len-2.0*h);
      if (y<y_1 && y>y_0){
        return true;
      }
      else{
        return false;
      }
    }
  }
  else if (shapeflag==1){
    y_1 = h*(zeta-1)*x/(2.0*len) +h/2.0;
    y_0 = -y_1;
    if (x>(-len) && x<(len) && y<y_1 && y>y_0){
      return true;
    }
    else{
      return false;
    }
  }
  else {
    float A = 1.0*0.2;
    y_1 = A/(2.0*len) - 2.0*A*x/(4.0*len*len);
    y_0 = -(A/(2.0*len) - 2.0*A*x/(4.0*len*len));
    if (y<y_1 && y>y_0){
      return true;
    }
    else{
      return false;
    }
  }


}
