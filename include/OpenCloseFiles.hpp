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
#include <stdlib.h>
#include <string>
#include <vector>

void OpenCloseFiles(char* buffer,std::vector<std::ofstream>& files,const int& realtimeplot,const int& endflag);


#endif
