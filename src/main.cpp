//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>
#include <vector>
//#include <stdlib.h>
#include "element.hpp"
#include "BoundToGrid1.hpp"
#include "BoundToGrid2.hpp"
#include "smallfunctions.hpp"
#include "GridToBound.hpp"
#include "NavierStokes.hpp"
#include "GlobalToLocal.hpp"
#include "OpenCloseFiles.hpp"
#include "OutputData.hpp"
#include "MatrixAdhesion.hpp"
#include "LocalToGlobal.hpp"
#include "AdjacentForces.hpp"
#include "ReadParameters.hpp"
#include "Initialise.hpp"

using namespace std;
using namespace arma;

// System parameters
int   Ng       ;         // Fluid grid size
int   Nb       ;         // Number of boundary points
float dims=1   ;         // Fluid grid dimensions
float cen =0   ;         // Fluid centre point
float Src =0   ;         // Source strength (= cell growth rate)
float rho      ;         // Fluid density
float mu       ;         // Fluid viscosity
float dt       ;         // Time step in seconds
float t_max    ;         // Max run time in seconds
float t_output ;         // Output interval in seconds
float tension  ;         // Cell cortical tension
float adhesion ;
float re       ;
float t     = 0;         // Run time in seconds
int   nloop = 0;         // Just counts how many time steps there have been so far
int   realtimeplot;      // Flag for real time plotting
int   plotfluid;
vector<ofstream> files;  // Set of output files
cube  xg       ;
mat   sg       ;
cube  fg       ;
cube  vg       ;
cube  ug       ;
mat   xbglobal ;
mat   ubglobal ;
mat   fbglobal ;
int   Nbs      = 4;
float xmin;
float xmax;
mat sb       = mat(2,Nbs,fill::zeros);
mat sbb      = mat(2,Nbs,fill::zeros);
vector<element> Elements;
float hg;
char outputfolder[26];

int main() {

  // Set up data output files
  OpenCloseFiles(outputfolder,files,realtimeplot,0);

  Initialise(files,Elements,Nb,Ng,rho,mu,re,tension,adhesion,dt,t_max,t_output,realtimeplot,plotfluid,xmin,xmax,hg,xg,sg,fg,vg,ug,xbglobal,ubglobal,fbglobal);

  OutputData(outputfolder,files,Elements,xbglobal,xg,fg,Nb,Ng,nloop,realtimeplot,1,plotfluid,xmin,xmax);

  // Iterate system over time
  cout << "Initialising shape evolution" << endl;
  while (t<t_max) {
    // Calculate tension forces
    AdjacentForces(Elements,hg,re,Nb,tension);
    // Calculate adhesion forces
    MatrixAdhesion(Elements,Nb,adhesion,hg,re);
    // Convert from local data to global arrays
    LocalToGlobal(Elements,xbglobal,fbglobal,ubglobal,Nb);
    // Calculate contributions of fluid sources
    BoundToGrid1(sg,Nbs,sb,sbb,xmin,xmax,hg,Ng);
    // Calculate contributions of boundary forces on fluid
    BoundToGrid2(fg,Nb,fbglobal,xbglobal,xmin,xmax,hg,Ng);
    // Compute grid velocity from NavierStokes
    NavierStokes(ug,vg,fg,sg,Ng,mu,rho,hg,dt);
    // Calculate contributions of fluid velocity to boundary point velocities
    GridToBound(ubglobal,xbglobal,vg,hg,xmin,xmax,Nb,Ng);
    // Update positions of boundary points according to calculated velocities
    xbglobal = xbglobal + dt*ubglobal;
    // Convert global arrays back to local data
    GlobalToLocal(Elements,xbglobal,Nb);
    // Increment time
    t = t+dt;
    // Write data to file at every output interval
    if (fmod(t,t_output)<dt){
      OutputData(outputfolder,files,Elements,xbglobal,xg,fg,Nb,Ng,nloop,realtimeplot,0,plotfluid,xmin,xmax);
      system("clear");
      printf("IBM progress: %f/%f\n",t,t_max);
    }
  }

  // Close data files
  OpenCloseFiles(outputfolder,files,realtimeplot,1);
  return 0;
}
