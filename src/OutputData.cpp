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
#include <armadillo>
#include "element.hpp"
#include <string>

using namespace std;
using namespace arma;

void OutputData(const char* buffer2,vector<ofstream>& files,const vector<element>& Elements,const mat& xbglobal,const cube& xg,const cube& fg,const int& Nb, const int& Ng,int& nloop,const int& realtimeplot,const int& startflag,const int& plotfluid,const float& xmin,const float& xmax){

  char  buffer[100];       // Dummy string for system calls

  for (int ii=0;ii<Nb;ii++){
    files[1] << Elements[ii].pos(0) << ", ";
    files[1] << Elements[ii].pos(1) << ", ";
    files[1] << Elements[ii].accumulatedEffector << ", ";
    files[1] << Elements[ii].n0 << ", ";
    files[1] << Elements[ii].n1 << endl;
  }
  files[2] << Nb << endl;
  files[1].flush();
  files[2].flush();

  if (startflag==1){
    // Write grid positions to file
    for (int ii=0;ii<Ng+1;ii++){
      files[3] << xg.slice(0).row(ii);
      files[4] << xg.slice(1).row(ii);
    }
    files[3].close();
    files[4].close();
  }

  if (plotfluid == 1){
    for (int ii=0;ii<Ng+1;ii++){
      files[5] << fg.slice(0).row(ii);
      files[6] << fg.slice(1).row(ii);
    }
    files[5].flush();
    files[6].flush();
  }

  if (realtimeplot==1) {
    // Call plotter
    sprintf(buffer,"python3 scripts/plottersingle.py %d %d %d %d %.3f %.3f ",nloop,Nb,Ng,2,xmin,xmax);
    system((string(buffer)+string(buffer2)).c_str());//system((string(buffer)+string(buffer2)+" &").c_str());
    nloop = nloop+1;
  }

}
