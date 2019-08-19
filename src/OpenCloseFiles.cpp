//
//  OpenFiles.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 08/02/2019.
//
//

#include "OpenCloseFiles.hpp"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>

using namespace std;

void OpenCloseFiles(char* buffer,vector<ofstream>& files,const int& realtimeplot,const int& endflag){

  vector<string> names;
  std::time_t result = std::time(nullptr);
  struct tm * timeinfo;
  timeinfo = localtime (&result);

  if (endflag  == 1){
    for (int ii=0;ii<files.size();ii++){
      files[ii].close();
    }
  }
  else{
    strftime (buffer,26, "output/%F-%H-%M-%S",timeinfo);
    system(("mkdir "+string(buffer)).c_str());
    names.push_back((string(buffer)+"/initialconditions.txt").c_str());
    names.push_back((string(buffer)+"/boundarypositions.txt").c_str());
    names.push_back((string(buffer)+"/nbounds.txt").c_str());
    names.push_back((string(buffer)+"/fluidvelocities0.txt").c_str());
    names.push_back((string(buffer)+"/fluidvelocities1.txt").c_str());
    names.push_back((string(buffer)+"/gridpositions0.txt").c_str());
    names.push_back((string(buffer)+"/gridpositions1.txt").c_str());
    files.resize(names.size());
    for (int ii=0;ii<names.size();ii++){
      files[ii].open(names[ii],ofstream::out);
    }
  }

}
