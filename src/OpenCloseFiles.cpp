//
//  OpenFiles.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 08/02/2019.
//
//

#include "OpenCloseFiles.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <stdio.h>

using namespace std;

void OpenCloseFiles(char*  buffer,vector<ofstream>& files,const int& endflag,const float& zeta,const float& adhesion,const float& conc){

  vector<string> names;
  std::time_t result = std::time(nullptr);
  struct tm * timeinfo;
  timeinfo = localtime (&result);


    sprintf (buffer, "output/%04d-%04d-%04d",static_cast<int>(zeta*1000),static_cast<int>(adhesion*1000),static_cast<int>(conc*1000));
    system(("mkdir "+string(buffer)).c_str());
    names.push_back((string(buffer)+"/initialconditions.txt").c_str());
    names.push_back((string(buffer)+"/boundarypositions.txt").c_str());
    names.push_back((string(buffer)+"/nbounds.txt").c_str());
    names.push_back((string(buffer)+"/gridpositions0.txt").c_str());
    names.push_back((string(buffer)+"/gridpositions1.txt").c_str());
    names.push_back((string(buffer)+"/fluidvelocities0.txt").c_str());
    names.push_back((string(buffer)+"/fluidvelocities1.txt").c_str());
    files.resize(names.size());
    for (int ii=0;ii<names.size();ii++){
      files[ii].open(names[ii],fstream::out);
    }


}

void OpenCloseFiles(char*  buffer,vector<ofstream>& files,const int& endflag,const int& realtimeplot){

  vector<string> names;
  std::time_t result = std::time(nullptr);
  struct tm * timeinfo;
  timeinfo = localtime (&result);

  if (endflag==1){
    for (int ii=0;ii<files.size();ii++){
      files[ii].close();
    }
    if (realtimeplot==1){
      system(("convert -delay 10 -loop 0 "+string(buffer)+"/plot*.png "+string(buffer)+"/plotanimated.gif").c_str());
    }
  }
}
