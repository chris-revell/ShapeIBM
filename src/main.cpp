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

using namespace std;
using namespace arma;

// System parameters
int   Numg        = 512;  // Fluid grid size
int   Nb          = 512;   // Number of boundary points
int   dims        = 10;   // Fluid grid dimensions
float cen         = 0;    // Fluid centre point
float Src         = 0.0;  // Source strength
float rho         = 1;    // Fluid density
float mu          = 10;   // Fluid viscosity
float len         = 2;    // Initial cell radius in micrometres
int   Numcells    = 1;    // Number of cells
float dt          = 1;   // Time step in seconds
float t           = 0;    // Run time in seconds
float t_max       = 200;  // Max run time in seconds
float t_output    = 10.0;  // Output interval in seconds
float tension     = 1;  // Cell cortical tension
int   nloop       = 0;    // Just counts how many time steps there have been so far
int   realtimeplot= 0;    // Flag for real time plotting
int   exitval;            // Dummy variable for system calls
char  buffer[50];         // Dummy string for system calls
vector<ofstream> files;   // Set of output files


int main() {
  //ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max,tension);
  exitval = system("rm output/grid*txt; rm output/fluid*txt;rm output/nbounds.txt;rm output/boundarypositions.txt;rm output/volume.txt;");//rm output/montageanimated.gif;");
  // Create whole tissue system
  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu,dt);
  // Add cell objects to tissue system
  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,0,0,tension,adhesion);
  }
  // Set up data output files
  OpenCloseFiles(files,realtimeplot,Tissue);
  // Iterate system over time
  while (t<t_max) {
    //Tissue.BoundaryRefinement();
    Tissue.UpdateSources();
    //-- Tension --
    for (int ii=0;ii<Tissue.Nc;ii++){
      AdjacentForces(Tissue.Cells[ii]);
    }
    //-- Adhesion
    for (int ii=0;ii<Tissue.Nc;ii++){
      MatrixAdhesion(Tissue.Cells[ii]);
    }
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
    GlobalToLocal(Tissue);
    // Write data to file
    if (fmod(t,t_output)<Tissue.dt){
      OutputData(files,t,Tissue,nloop,realtimeplot);
      printf("%f/%f\n",t,t_max);
    }
    // Increment time
    t = t+Tissue.dt;
  }
  // Close data files
  OpenCloseFiles(files,realtimeplot,Tissue);
  return 0;
}
