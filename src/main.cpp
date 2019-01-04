//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <vector>
#include <stdlib.h>
#include "cell.hpp"
#include "element.hpp"
#include "tissue.hpp"
#include "BoundToGrid1.hpp"
#include "BoundToGrid2.hpp"
#include "smallfunctions.hpp"
#include "GridToBound.hpp"
#include "NavierStokes.hpp"
#include "GlobalToLocal.hpp"
#include "ReadParams.hpp"

using namespace std;
using namespace arma;

// System parameters
int   Numg     = 512;   // Fluid grid size
int   Nb       = 64;    // Number of boundary points
int   dims     = 10;    // Fluid grid dimensions
float cen      = 0;     // Fluid centre point
float Src      = 0.0;   // Source strength
float rho      = 1;     // Fluid density
float mu       = 5;     // Fluid viscosity
float xi       = 0.01;  // Stochastic magnitude
float len      = 1;     // Initial cell radius in micrometres
int   Numcells = 1;     // Number of cells
float dt       = 1;     // Time step in seconds
float t        = 0;     // Run time in seconds
float t_max    = 2000;  // Max run time in seconds
float t_output = 10.0;  // Output interval in seconds
float tension  = 0.1;   // Cell cortical tension
int   nloop    = 0;     // Just counts how many time steps there have been so far
int   exitval;          // Dummy variable for system calls
char  buffer[50];       // Dummy string for system calls


int main() {

  //ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max,tension);
  exitval = system("rm output/velocity*.png;rm output/velocityanimated.gif;");

  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu,xi,dt);

  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,0,0,tension);
  }
  Tissue.UpdateSources();
  Tissue.CombineBoundaries();

  ofstream file1;
  file1.open ("output/boundarypositions.txt", ofstream::out);
  ofstream file2;
  file2.open ("output/nbounds.txt", ofstream::out);
  ofstream file3;
  file3.open ("output/volume.txt", ofstream::out);
  ofstream file4;
  file4.open ("output/fluidvelocities0.txt", ios::out);
  ofstream file5;
  file5.open ("output/fluidvelocities1.txt", ios::out);
  ofstream file6;
  file6.open ("output/gridpositions0.txt", ios::out);
  ofstream file7;
  file7.open ("output/gridpositions1.txt", ios::out);
  // Write grid positions to file
  for (int ii=0;ii<Numg+1;ii++){
    file6 << Tissue.xg.slice(0).row(ii);
    file7 << Tissue.xg.slice(1).row(ii);
  }
  // Write initial data to file //
  for (int ii=0;ii<Tissue.Nb;ii++){
    file1 << Tissue.xbglobal(0,ii) << ", ";
    file1 << Tissue.xbglobal(1,ii) << endl;
  }
  file1.flush();

  while (t<t_max) {

    //Tissue.BoundaryRefinement();
    Tissue.UpdateSources();
    //-- boundary forces --//
    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].AdjacentForces();
    }
    Tissue.CombineBoundaries();
    Tissue.ubglobal.zeros();
    //-- grid sources --//
    BoundToGrid1(Tissue);
    //-- grid forces --//
    BoundToGrid2(Tissue);
    //-- compute grid velocity from NavierStokes --//
    NavierStokes(Tissue);
    Tissue.ug = Tissue.vg;
    //-- boundary velocities --//
    GridToBound(Tissue);
    //-- new position of boundary points --//
    Tissue.UpdatePositions();
    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].UpdateCom();
    }
    GlobalToLocal(Tissue);

    if (fmod(t,t_output)<Tissue.dt){
      // Write data to file //
      cout << t << endl;
      for (int ii=0;ii<Tissue.Nb;ii++){
        file1 << Tissue.xbglobal(0,ii) << ", ";
        file1 << Tissue.xbglobal(1,ii) << endl;
      }
      file2 << Tissue.Nb << endl;
      file3 << t << " " << Tissue.Cells[0].CalculateVolume() << endl;
      for (int ii=0;ii<Numg+1;ii++){
        file4 << Tissue.vg.slice(0).row(ii);
        file5 << Tissue.vg.slice(1).row(ii);
      }
      file1.flush();
      file2.flush();
      file3.flush();
      file4.flush();
      file5.flush();
      // Call plotter
      exitval = sprintf(buffer,"python3 scripts/velocityplottersingle.py %d %d %d %d &",nloop,Tissue.Nb,Tissue.Ng,1);
      exitval = system(buffer);
      nloop = nloop+1;
    }
    //printf("%f/%f\n",t,t_max);

    t = t+Tissue.dt;

  }   // for loop_num
  file1.close();
  file2.close();
  file3.close();
  file4.close();
  file5.close();
  exitval = system("convert -delay 10 -loop 0 output/velocitytest*.png output/velocityanimated.gif");
  return 0;
}
