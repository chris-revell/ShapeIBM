//
//  ReadParameters.cpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 20/03/2019.
//
//

#include "ReadParameters.hpp"

//
//  ReadParameters.cpp
//  anisotropic-division
//
//  Created by Christopher Revell on 15/03/2019.
//
//

#include "ReadParameters.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void ReadParameters(int& Numg,int& Nb,float& dims,float& cen,float& Src,float& rho,float& mu,float& len,int& Numcells,float& dt,float& t_max,float& t_output,float& tension,float& adhesion,float& diffusionconstant,int& realtimeplot){
	static const std::streamsize max = std::numeric_limits<std::streamsize>::max();
	std::vector<float> values;
	string input;
	float value;
	ifstream infile("input/parameters.txt");

	while(infile.ignore(max,' ') >> input >> value)
	{
		values.push_back(value);
	}

  Numg               = values[0];
  Nb                 = values[1];
  dims               = values[2];
  cen                = values[3];
  Src                = values[4];
  rho                = values[5];
  mu                 = values[6];
  len                = values[7];
  Numcells           = values[8];
  dt                 = values[9];
  t_max              = values[10];
  t_output           = values[11];
  tension            = values[12];
  adhesion           = values[13];
  diffusionconstant  = values[14];
  realtimeplot       = values[15];

	cout << "Numg               " << Numg               << endl;
  cout << "Nb                 " << Nb                 << endl;
  cout << "dims               " << dims               << endl;
  cout << "cen                " << cen                << endl;
  cout << "Src                " << Src                << endl;
  cout << "rho                " << rho                << endl;
  cout << "mu                 " << mu                 << endl;
  cout << "len                " << len                << endl;
  cout << "Numcells           " << Numcells           << endl;
  cout << "dt                 " << dt                 << endl;
  cout << "t_max              " << t_max              << endl;
  cout << "t_output           " << t_output           << endl;
  cout << "tension            " << tension            << endl;
  cout << "adhesion           " << adhesion           << endl;
  cout << "diffusionconstant  " << diffusionconstant  << endl;
  cout << "realtimeplot       " << realtimeplot       << endl;

}
