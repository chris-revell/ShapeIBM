//
//  OpenFiles.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 08/02/2019.
//
//

#ifndef OpenFiles_hpp
#define OpenFiles_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>

void OpenCloseFiles(char*  buffer,std::vector<std::ofstream>& files,const int& endflag,const float& zeta,const float& adhesion,const float& conc,const float& len);
void OpenCloseFiles(char*  buffer,std::vector<std::ofstream>& files,const int& endflag,const int& realtimeplot);

#endif
