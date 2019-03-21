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
//#include "ReadParams.hpp"
#include "OpenCloseFiles.hpp"
#include "OutputData.hpp"
#include "MatrixAdhesion.hpp"
#include "LocalToGlobal.hpp"
#include "AdjacentForces.hpp"
#include <ReadParameters.hpp>

using namespace std;
using namespace arma;

// System parameters
int   Numg              = 256;   // Fluid grid size
int   Nb                = 256;   // Number of boundary points
float dims              = 10;    // Fluid grid dimensions
float cen               = 10;    // Fluid centre point
float Src               = 0.0;   // Source strength (= cell growth rate)
float rho               = 1;     // Fluid density
float mu                = 10;    // Fluid viscosity
float len               = 2;     // Typical cell radius in micrometres
int   Numcells          = 1;     // Number of cells
float dt                = 1;     // Time step in seconds

float t_max             = 80000; // Max run time in seconds
float t_output          = 500.0; // Output interval in seconds
float tension           = 1.0;   // Cell cortical tension
float adhesion          = 0.5;
float diffusionconstant = 0.03;

float t                 = 0;     // Run time in seconds
int   nloop       = 0;     // Just counts how many time steps there have been so far
// System control flags
int   realtimeplot= 1;     // Flag for real time plotting
//int   exitval;           // Dummy variable for system calls
char  buffer[50];          // Dummy string for system calls
vector<ofstream> files;    // Set of output files


int main() {
  //ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max,tension);
  system("rm output/grid*txt; rm output/fluid*txt;rm output/nbounds.txt;rm output/boundarypositions.txt;rm output/volume.txt;");//rm output/montageanimated.gif;");

  ReadParameters(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,dt,t_max,t_output,tension,adhesion,diffusionconstant,realtimeplot);

  // Create whole tissue system
  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu,dt,cen);
  // Add cell objects to tissue system
  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,10,10,tension,adhesion);
  }
  // Set up data output files
  OpenCloseFiles(files,realtimeplot,Tissue);
  // Iterate system over time
  while (t<t_max) {
    //Tissue.BoundaryRefinement();
    Tissue.UpdateSources();
    //-- Tension --
    //Tissue.Cells[0].ctension = t*tension/t_max;
    for (int ii=0;ii<Tissue.Nc;ii++){
      AdjacentForces(Tissue.Cells[ii],t,diffusionconstant);
    }
    //-- Adhesion
    MatrixAdhesion(Tissue,diffusionconstant,t);
    //-- Local and global array book keeping
    LocalToGlobal(Tissue);
    //-- grid sources
    BoundToGrid1(Tissue);
    //-- grid forces
    BoundToGrid2(Tissue);
    //-- compute grid velocity from NavierStokes
    NavierStokes(Tissue);
    //-- boundary velocities
    GridToBound(Tissue);
    //-- new position of boundary points
    Tissue.UpdatePositions();
    //-- Update cell centres of mass
    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].UpdateCom();
    }
    //-- Local and global array book keeping
    GlobalToLocal(Tissue);
    // Write data to file
    if (fmod(t,t_output)<Tissue.dt){
      OutputData(files,t,Tissue,nloop,realtimeplot,diffusionconstant);
      printf("%f/%f\n",t,t_max);
    }
    // Increment time
    t = t+Tissue.dt;
  }
  // Write data to file
  if (fmod(t,t_output)<Tissue.dt){
    OutputData(files,t,Tissue,nloop,realtimeplot,diffusionconstant);
    printf("%f/%f\n",t,t_max);
  }
  // Close data files
  OpenCloseFiles(files,realtimeplot,Tissue);
  return 0;
}
