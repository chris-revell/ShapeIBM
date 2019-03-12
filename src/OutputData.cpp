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

void OutputData(vector<ofstream>& files,const float& t,tissue& Tissue,int& nloop,const int& realtimeplot,const float& diffusionconstant)  {
  //int   exitval;          // Dummy variable for system calls
  char  buffer[75];       // Dummy string for system calls
  for (int ii=0;ii<Tissue.Nb;ii++){
    files[0] << Tissue.xbglobal(0,ii) << ", ";
    files[0] << Tissue.xbglobal(1,ii) << endl;//", ";
    //files[0] << Tissue.fbglobal(0,ii) << ", ";
    //files[0] << Tissue.fbglobal(1,ii) << endl;
  }
  files[1] << Tissue.Nb << endl;
  files[2] << t << " " << Tissue.Cells[0].CalculateVolume() << endl;
  for (int ii=0;ii<Tissue.Ng+1;ii++){
    files[3] << Tissue.fg.slice(0).row(ii);
    files[4] << Tissue.fg.slice(1).row(ii);
  }
  files[0].flush();
  files[1].flush();
  files[2].flush();
  files[3].flush();
  files[4].flush();
  if (realtimeplot==1) {
    // Call plotter
    sprintf(buffer,"python3 scripts/plottersingle.py %d %d %d %d %.4f %.4f &",nloop,Tissue.Nb,Tissue.Ng,2,t,diffusionconstant);
    system(buffer);
    nloop = nloop+1;
  }
}
