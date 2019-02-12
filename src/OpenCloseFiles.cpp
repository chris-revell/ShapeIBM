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

using namespace std;

void OpenCloseFiles(vector<ofstream>& files,int& realtimeplot){
  int   exitval;          // Dummy variable for system calls
  vector<string> names;// = {"boundarypositions.txt","nbounds.txt","volume.txt","fluidvelocities0.txt","fluidvelocities1.txt","gridpositions0.txt","gridpositions1.txt"};

  names.push_back("boundarypositions.txt");
  names.push_back("nbounds.txt");
  names.push_back("volume.txt");
  names.push_back("fluidvelocities0.txt");
  names.push_back("fluidvelocities1.txt");
  names.push_back("gridpositions0.txt");
  names.push_back("gridpositions1.txt");

  files.resize(names.size());
  if (files[0].is_open()){
    for (int ii=0;ii<files.size();ii++){
      files[ii].close();
    }
    if (realtimeplot==1){
      exitval = system("convert -delay 10 -loop 0 output/velocitytest*.png output/velocityanimated.gif");
    }
  }
  else{
    for (int ii=0;ii<names.size();ii++){
      files[ii].open("output/"+names[ii],ofstream::out);
    }
  }

}
