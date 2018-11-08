//
//  ReadParams.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 01/11/2018.
//
//

#include "ReadParams.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>

using namespace std;

void ReadParams(int& Numg,int& Nb,int& dims,float& cen,float& Src,float& rho,float& mu,float& dt,float& len,int& NumLoop,int& Numcells){
  ifstream infile("input/parameters.txt");
  string line;
  vector<string> params;

  if (infile.is_open()) {

    while (std::getline(infile,line)){
      params.push_back(line);
    }
    Numg  = stoi(params[0]);
    Nb = stoi(params[1]);
    dims = stoi(params[2]);
    cen = stof(params[3]);
    Src = stof(params[4]);
    rho = stof(params[5]);
    mu = stof(params[6]);
    dt = stof(params[7]);
    len = stof(params[8]);
    NumLoop = stoi(params[9]);
    Numcells = stoi(params[10]);
  }else{
    cout << "Parameter file input/parameters.txt not found" << endl;
  }
}