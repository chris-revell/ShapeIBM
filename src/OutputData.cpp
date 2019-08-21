//
//  OutputData.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 10/02/2019.
//
//

#include "OutputData.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <armadillo>

using namespace std;
using namespace arma;

void OutputData(const char* buffer2,vector<ofstream>& files,const vector<element>& Elements,const mat& xbglobal,const cube& xg,const cube& fg,const float& t,const int& Nb, const int& Ng,int& nloop,const int& realtimeplot,const int& startflag,const float& xmin,const float& xmax){

  char  buffer[100];       // Dummy string for system calls

  for (int ii=0;ii<Nb;ii++){
    files[1] << Elements[ii].pos(0) << ", ";
    files[1] << Elements[ii].pos(1) << ", ";
    files[1] << Elements[ii].accumulatedEffector << endl;

    //files[1] << Tissue.fbglobal(0,ii) << ", ";
    //files[1] << Tissue.fbglobal(1,ii) << endl;
  }
  files[2] << Nb << endl;
  for (int ii=0;ii<Ng+1;ii++){
    files[3] << fg.slice(0).row(ii);
    files[4] << fg.slice(1).row(ii);
  }
  files[1].flush();
  files[2].flush();
  files[3].flush();
  files[4].flush();
  if (startflag<1 && startflag>=0){
    // Write grid positions to file
    for (int ii=0;ii<Ng+1;ii++){
      files[5] << xg.slice(0).row(ii);
      files[6] << xg.slice(1).row(ii);
    }
    files[5].close();
    files[6].close();
  }
  if (realtimeplot==1) {
    // Call plotter
    sprintf(buffer,"python3 scripts/plottersingle.py %d %d %d %d %.3f %.3f %.4f ",nloop,Nb,Ng,2,xmin,xmax,t);
    system((string(buffer)+string(buffer2)).c_str());//system((string(buffer)+string(buffer2)+" &").c_str());
    nloop = nloop+1;
  }
  if (startflag<0 && realtimeplot==1){
    system(("convert -delay 10 -loop 0 "+string(buffer2)+"/plot*.png "+string(buffer2)+"/plotanimated.gif &").c_str());
  }
}
