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
#include <tissue.hpp>

using namespace std;

void OutputData(char* buffer2,vector<ofstream>& files,const float& t,tissue& Tissue,int& nloop,const int& realtimeplot,const int& startflag)  {
  //int   exitval;        // Dummy variable for system calls
  char  buffer[100];       // Dummy string for system calls
  for (int ii=0;ii<Tissue.Nb;ii++){
    files[1] << Tissue.xbglobal(0,ii) << ", ";
    files[1] << Tissue.xbglobal(1,ii) << ", ";
    files[1] << Tissue.Elements[ii].accumulatedEffector << endl;
    //files[1] << Tissue.fbglobal(0,ii) << ", ";
    //files[1] << Tissue.fbglobal(1,ii) << endl;
  }
  files[2] << Tissue.Nb << endl;
  for (int ii=0;ii<Tissue.Ng+1;ii++){
    files[3] << Tissue.fg.slice(0).row(ii);
    files[4] << Tissue.fg.slice(1).row(ii);
  }
  files[1].flush();
  files[2].flush();
  files[3].flush();
  files[4].flush();
  if (startflag<1 && startflag>=0){
    // Write grid positions to file
    for (int ii=0;ii<Tissue.Ng+1;ii++){
      files[5] << Tissue.xg.slice(0).row(ii);
      files[6] << Tissue.xg.slice(1).row(ii);
    }
    files[5].close();
    files[6].close();
  }
  if (realtimeplot==1) {
    // Call plotter
    sprintf(buffer,"python3 scripts/plottersingle.py %d %d %d %d %.3f %.3f %.4f ",nloop,Tissue.Nb,Tissue.Ng,2,Tissue.xmin,Tissue.xmax,t);    
    system((string(buffer)+string(buffer2)+" &").c_str());
    nloop = nloop+1;
  }
  if (startflag<0 && realtimeplot==1){
    system(("convert -delay 10 -loop 0 "+string(buffer2)+"/plot*.png "+string(buffer2)+"/plotanimated.gif &").c_str());
  }
}
