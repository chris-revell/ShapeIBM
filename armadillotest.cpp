#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void arraytest(mat& testarray,const int& x, const int&y){
  cout << testarray(x,y) << endl;
  testarray(x,y) = 1;
}

//void updatearray(cx_cube& data,const int& x, const int&y,const int& z,const float& dx,const float& dy,const float& dz){
//  data(x,y,z) = data(x,y,z)+ dx*dy*dz;
//}

int main() {
  const int Ng = 64;
  mat fvgg     = mat(Ng+1,Ng+1,fill::zeros);
  cx_cube data     = cx_cube(7,8,9,fill::randu);
  cx_mat vg1   = cx_mat(Ng,Ng,fill::zeros);

  //cout << fvgg(1,1) << endl;
  //cout << vg1(1,1) << endl;
  //fvgg(1,1) = 3.1415;
  //cout << fvgg(1,1) << endl;
  //cout << size(fvgg) << endl;
  //cout << size(fvgg(span(0,Ng-1),span(0,Ng-1))) << endl;
  //vg1 = fft2(fvgg(span(0,Ng-1),span(0,Ng-1)));
  //cout << vg1(1,1) << endl;

  //fvgg(span(0,Ng-1),span(0,Ng-1)) = real(vg1);
  //cout << fvgg(1,1) << endl;
  //cout << fvgg(Ng,Ng) << endl;

  arraytest(fvgg,5,10);
  //cout << fvgg(5,10) << endl;

  //cout << data(5,6,7) << endl;
  ////updatearray(data,5,6,7,1,2,3);
  //cout << data(5,6,7) << endl;

  data(5,1,1) = 1;

  //cout << size(data.slice(1)) << endl;
  //cout << data.slice(1);


  cout << "#######" << endl;
  cout << data(1,2,3) << endl;
  cout << real(data(1,2,3)) << endl;

}
