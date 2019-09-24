//
//  ReadParameters.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 20/03/2019.
//
//

#include "ReadParameters.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void ReadParameters(int argc,char *argv[],ofstream& file,int& Ng,float& rho,float& mu,float& len,float& h,float& zeta,float& re,float& tension,float& adhesion,float& D,float& conc,float& tdif_max,float& dt,float& t_max,float& t_output,int& realtimeplot,int& plotfluid,int& shapeflag){
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


	file << "Ng                 " << Ng              << endl;
  file << "rho                " << rho             << endl;
  file << "mu                 " << mu              << endl;
  file << "len                " << len             << endl;
  file << "h                  " << h               << endl;
	file << "zeta               " << zeta            << endl;
  file << "re                 " << re              << endl;
  file << "tension            " << tension         << endl;
  file << "adhesion           " << adhesion        << endl;
  file << "D                  " << D               << endl;
  file << "conc               " << conc             << endl;
  file << "tdif_max           " << tdif_max        << endl;
  file << "dt                 " << dt              << endl;
  file << "t_max              " << t_max           << endl;
  file << "t_output           " << t_output        << endl;
  file << "realtimeplot       " << realtimeplot    << endl;
	file << "plotfluid          " << plotfluid       << endl;
  file << "shapeflag          " << shapeflag       << endl;
	file.flush();

}
