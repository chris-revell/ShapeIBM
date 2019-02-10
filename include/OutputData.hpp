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
#include <tissue.hpp>

void OutputData(std::vector<std::ofstream>& files,float& t,tissue& Tissue);

#endif
