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
//#include <mkl.h>
//#include <omp.h>
#include "OpenCloseFiles.hpp"
#include "OutputData.hpp"

using namespace std;
using namespace arma;

// System parameters
int   Numg     = 64;   // Fluid grid size
int   Nb       = 32;    // Number of boundary points
int   dims     = 10;    // Fluid grid dimensions
float cen      = 0;     // Fluid centre point
float Src      = 0.0;   // Source strength
float rho      = 1;     // Fluid density
float mu       = 10;     // Fluid viscosity
//float xi       = 0.0;  // Stochastic magnitude
float len      = 2;     // Initial cell radius in micrometres
int   Numcells = 1;     // Number of cells
float dt       = 1;     // Time step in seconds
float t        = 0;     // Run time in seconds
float t_max    = 200;  // Max run time in seconds
float t_output = 1.0;  // Output interval in seconds
float tension  = 0.1;   // Cell cortical tension
int   nloop    = 0;     // Just counts how many time steps there have been so far
int   exitval;          // Dummy variable for system calls
char  buffer[50];       // Dummy string for system calls
vector<ofstream> files;


int main() {

  //ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max,tension);
  //exitval = system("rm output/velocity*.png;rm output/velocityanimated.gif;rm grid*txt; rm fluid*txt");
  //exitval = system("rm output/montage*.png;rm output/volume*.png;rm output/montageanimated.gif;rm grid*txt; rm fluid*txt");

  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu,dt);

  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,0,0,tension);
  }

  OpenCloseFiles(files);

  // Write grid positions to file
  for (int ii=0;ii<Numg+1;ii++){
    files[5] << Tissue.xg.slice(0).row(ii);
    files[6] << Tissue.xg.slice(1).row(ii);
  }

  while (t<t_max) {

    //Tissue.BoundaryRefinement();
    Tissue.UpdateSources();
    //-- boundary forces --//
    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].AdjacentForces();
    }
    Tissue.CombineBoundaries();
    if (t<100){
      Tissue.MatrixAdhesions();
    }
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
      OutputData(files,t,Tissue);
      // Call plotter
      //exitval = sprintf(buffer,"python3 scripts/velocityplottersingle.py %d %d %d %d &",nloop,Tissue.Nb,Tissue.Ng,1);
      //exitval = system(buffer);
      //nloop = nloop+1;
    }
    printf("%f/%f\n",t,t_max);

    t = t+Tissue.dt;

  }   // for loop_num

  OpenCloseFiles(files);
  //exitval = system("convert -delay 10 -loop 0 output/velocitytest*.png output/velocityanimated.gif");
  return 0;
}
