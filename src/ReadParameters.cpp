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

void ReadParameters(ofstream& file,int& Ng,float& rho,float& mu,float& len,float& h,float& zeta,float& re,float& tension,float& adhesion,float& D,float& conc,float& tdif_max,float& dt,float& t_max,float& t_output,int& realtimeplot){
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

	cout << "Ng                 " << Ng              << endl;
  cout << "rho                " << rho             << endl;
  cout << "mu                 " << mu              << endl;
  cout << "len                " << len             << endl;
  cout << "h                  " << h               << endl;
	cout << "zeta               " << zeta            << endl;
  cout << "re                 " << re              << endl;
  cout << "tension            " << tension         << endl;
  cout << "adhesion           " << adhesion        << endl;
  cout << "D                  " << D               << endl;
  cout << "tdif_max           " << tdif_max        << endl;
  cout << "dt                 " << dt              << endl;
  cout << "t_max              " << t_max           << endl;
  cout << "t_output           " << t_output        << endl;
  cout << "realtimeplot       " << realtimeplot    << endl;

}
