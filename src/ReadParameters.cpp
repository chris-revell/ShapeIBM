//
//  ReadParameters.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 20/03/2019.
//
//

#include "ReadParameters.hpp"
#include "OpenCloseFiles.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void ReadParameters(int argc,char *argv[],char* outputfolder,vector<ofstream>& files,int& Ng,float& rho,float& mu,float& len,float& h,float& zeta,float& re,float& tension,float& adhesion,float& D,float& conc,float& tdif_max,float& dt,float& t_max,float& t_output,int& realtimeplot,int& plotfluid,int& shapeflag){
	//static const std::streamsize max = std::numeric_limits<std::streamsize>::max();
//	std::vector<float> values;
//	string input;
//	float value;
//
//	ifstream infile("input/parameters.txt");
//
//	while(infile.ignore(max,' ') >> input >> value)
//	{
//		values.push_back(value);
//	}
//
  Ng            = atoi(argv[1]);
  rho           = atof(argv[2]);
  mu            = atof(argv[3]);
  len           = atof(argv[4]);
  h             = atof(argv[5]);
	zeta					= atof(argv[6]);
  re            = atof(argv[7]);
  tension       = atof(argv[8]);
  adhesion      = atof(argv[9]);
  D             = atof(argv[10]);
	conc 					= atof(argv[11]);
  tdif_max      = atof(argv[12]);
  dt            = atof(argv[13]);
  t_max         = atof(argv[14]);
	t_output      = atof(argv[15]);
	realtimeplot  = atoi(argv[16]);
	plotfluid     = atoi(argv[17]);
  shapeflag     = atoi(argv[18]);

  OpenCloseFiles(outputfolder,files,0,zeta,adhesion,conc,len);


	files[0] << "Ng                 " << Ng              << endl;
  files[0] << "rho                " << rho             << endl;
  files[0] << "mu                 " << mu              << endl;
  files[0] << "len                " << len             << endl;
  files[0] << "h                  " << h               << endl;
	files[0] << "zeta               " << zeta            << endl;
  files[0] << "re                 " << re              << endl;
  files[0] << "tension            " << tension         << endl;
  files[0] << "adhesion           " << adhesion        << endl;
  files[0] << "D                  " << D               << endl;
  files[0] << "conc               " << conc            << endl;
  files[0] << "tdif_max           " << tdif_max        << endl;
  files[0] << "dt                 " << dt              << endl;
  files[0] << "t_max              " << t_max           << endl;
  files[0] << "t_output           " << t_output        << endl;
  files[0] << "realtimeplot       " << realtimeplot    << endl;
	files[0] << "plotfluid          " << plotfluid       << endl;
  files[0] << "shapeflag          " << shapeflag       << endl;
	files[0].flush();

}
