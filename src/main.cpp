//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>
#include <vector>
#include <stdlib.h>
//#include "cell.hpp"
#include "element.hpp"
#include "tissue.hpp"
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

int realtimeplot;
float t_max,t_output;
float t     = 0;     // Run time in seconds
int   nloop = 0;     // Just counts how many time steps there have been so far
//char  buffer[50];          // Dummy string for system calls
char outputfolder[26];
vector<ofstream> files;    // Set of output files

int main() {

  // Set up data output files
  OpenCloseFiles(outputfolder,files,realtimeplot,0);

  tissue Tissue = Initialise(realtimeplot,t_max,t_output,files[0]);

  // Iterate system over time
  cout << "Initialising shape evolution" << endl;
  while (t<t_max) {
    // Calculate tension forces
    AdjacentForces(Tissue);
    // Calculate adhesion forces
    MatrixAdhesion(Tissue);
    // Convert from local data to global arrays
    LocalToGlobal(Tissue);
    // Calculate contributions of fluid sources
    BoundToGrid1(Tissue);
    // Calculate contributions of boundary forces on fluid
    BoundToGrid2(Tissue);
    // Compute grid velocity from NavierStokes
    NavierStokes(Tissue);
    // Calculate contributions of fluid velocity to boundary point velocities
    GridToBound(Tissue);
    // Update positions of boundary points according to calculated velocities
    Tissue.UpdatePositions();
    // Convert global arrays back to local data
    GlobalToLocal(Tissue);
    // Write data to file at every output interval
    if (fmod(t,t_output)<Tissue.dt){
      OutputData(outputfolder,files,t,Tissue,nloop,realtimeplot,static_cast<int>(t+1-Tissue.dt));
    }
    // Increment time
    t = t+Tissue.dt;
    system("clear");
    printf("IBM progress: %f/%f\n",t,t_max);
  }

  // Write data to file at final system state
  OutputData(outputfolder,files,t,Tissue,nloop,realtimeplot,-1);
  printf("%f/%f\n",t,t_max);
  // Close data files
  OpenCloseFiles(outputfolder,files,realtimeplot,1);
  return 0;
}
