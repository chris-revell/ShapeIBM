//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <vector>
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
int   Numg=512;            // fluid grid size
int   Nb  =256;            // number of boundary points
int   dims=10;             // Fluid grid dimensions
float cen=0;              // Fluid centre point
float Src=0.0;              // source strength
float rho=1;              // fluid density
float mu=1;               // fluid viscosity
float len=1;            // Initial cell radius in micrometres
int   Numcells=1;         // number of cells
float t=0;                // Run time in seconds
float t_max=10;            // Max run time in seconds
float corticaltension=0.001;// Cell cortical tension

int main() {

  //ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max,corticaltension);

  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu);

  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,3*len+3*ii*len,0,corticaltension);
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

    if (fmod(t,0.1)<Tissue.dt){
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
    }
    //printf("%f/%f\n",t,t_max);

    t = t+Tissue.dt;

  }   // for loop_num
  for (int ii=0;ii<Numg+1;ii++){
    file6 << Tissue.xg.slice(0).row(ii);
    file7 << Tissue.xg.slice(1).row(ii);
  }
  file1.close();
  file2.close();
  file3.close();
  file4.close();
  file5.close();
  return 0;
}
