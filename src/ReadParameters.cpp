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

void ReadParameters(ofstream& file,int& Numg,int& Nb,float& dims,float& cen,float& Src,float& rho,float& mu,float& len,float& dt,float& t_max,float& t_output,float& tension,float& adhesion,int& realtimeplot,float& h,float& A,float& alpha,float& D,float& tdif_max,float& re){
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
  dt                 = values[8];
  t_max              = values[9];
  t_output           = values[10];
  tension            = values[11];
  adhesion           = values[12];
  realtimeplot       = values[13];
	h									 = values[14];
	A									 = values[15];
	alpha							 = values[16];
	D									 = values[17];
	tdif_max					 = values[18];
	re								 = values[19];

	file << "Numg               " << Numg               << endl;
  file << "Nb                 " << Nb                 << endl;
  file << "dims               " << dims               << endl;
  file << "cen                " << cen                << endl;
  file << "Src                " << Src                << endl;
  file << "rho                " << rho                << endl;
  file << "mu                 " << mu                 << endl;
  file << "len                " << len                << endl;
  file << "dt                 " << dt                 << endl;
  file << "t_max              " << t_max              << endl;
  file << "t_output           " << t_output           << endl;
  file << "tension            " << tension            << endl;
  file << "adhesion           " << adhesion           << endl;
  file << "realtimeplot       " << realtimeplot       << endl;
  file << "h                  " << h    			        << endl;
  file << "A                  " << A    			        << endl;
  file << "alpha              " << alpha			        << endl;
	file << "D			            " << D				          << endl;
	file << "tdif_max           " << tdif_max           << endl;
	file << "re			            " << re			            << endl;

}
