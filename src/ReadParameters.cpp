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

void ReadParameters(ofstream& file,int& Ng,float& rho,float& mu,float& len,float& h,float& zeta,float& re,float& tension,float& adhesion,float& D,float& conc,float& tdif_max,float& dt,float& t_max,float& t_output,int& realtimeplot,int& plotfluid){
	static const std::streamsize max = std::numeric_limits<std::streamsize>::max();
	std::vector<float> values;
	string input;
	float value;

	ifstream infile("input/parameters.txt");

	while(infile.ignore(max,' ') >> input >> value)
	{
		values.push_back(value);
	}

  Ng            = values[0];
  rho           = values[1];
  mu            = values[2];
  len           = values[3];
  h             = values[4];
	zeta					= values[5];
  re            = values[6];
  tension       = values[7];
  adhesion      = values[8];
  D             = values[9];
	conc 					= values[10];
  tdif_max      = values[11];
  dt            = values[12];
  t_max         = values[13];
	t_output      = values[14];
	realtimeplot  = values[15];
	plotfluid     = values[16];

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
  file << "tdif_max           " << tdif_max        << endl;
  file << "dt                 " << dt              << endl;
  file << "t_max              " << t_max           << endl;
  file << "t_output           " << t_output        << endl;
  file << "realtimeplot       " << realtimeplot    << endl;
	file << "plotfluid          " << plotfluid       << endl;
	file.flush();

}
