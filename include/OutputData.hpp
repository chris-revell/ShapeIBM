//
//  OutputData.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 10/02/2019.
//
//

#ifndef OutputData_hpp
#define OutputData_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <armadillo>
#include <element.hpp>

void OutputData(const char* buffer2,std::vector<std::ofstream>& files,const std::vector<element>& Elements,const arma::mat& xbglobal,const arma::cube& xg,const arma::cube& fg,const float& t,const int& Nb, const int& Ng,int& nloop,const int& realtimeplot,const int& startflag,const float& xmin,const float& xmax);

#endif
