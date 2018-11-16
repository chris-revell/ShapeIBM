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
int   Numg;    // fluid grid size
int   Nb;      // number of boundary points
int   dims;    // Fluid grid dimensions
float cen;     // Fluid centre point
float Src;     // source strength
float rho;     // fluid density
float mu;      // fluid viscosity
float len;     // Initial cell radius in micrometres
int   Numcells;// number of cells
float t;    // Run time in seconds
float t_max;// Max run time in seconds

int main() {

  t=0;
  ReadParams(Numg,Nb,dims,cen,Src,rho,mu,len,Numcells,t_max);

  tissue Tissue = tissue(Numg,dims,Nb,Src,rho,mu);

  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,3*ii*len,0);
  }
  Tissue.UpdateSources();
  Tissue.CombineBoundaries();

  ofstream file1;
  file1.open ("output/boundarypositions.txt", ofstream::out);
  ofstream file2;
  file2.open ("output/nbounds.txt", ofstream::out);
  // Write initial data to file //
  for (int ii=0;ii<Tissue.Nb;ii++){
    file1 << Tissue.xbglobal(0,ii) << ", ";
    file1 << Tissue.xbglobal(1,ii) << endl;
  }
  file1.flush();

  while (t<t_max) {

    Tissue.BoundaryRefinement();

    file2 << Tissue.Nb << endl;

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
    NavierStokes(Tissue.vg,Tissue.ug,Tissue.fg,Tissue.sg,Tissue.Ng,Tissue.rho,Tissue.mu,Tissue.dt,Tissue.hg);

    Tissue.ug = Tissue.vg;
    //-- boundary velocities --//
    GridToBound(Tissue.ubglobal,Tissue.xbglobal,Tissue.Nb,Tissue.vg,Tissue.Ng,Tissue.hg,Tissue.hg,Tissue.xmin,Tissue.xmax);

    //-- new position of boundary points --//
    Tissue.UpdatePositions();

    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].UpdateCom();
    }

    GlobalToLocal(Tissue);

    Tissue.UpdateSources();

    // Write data to file //
    for (int ii=0;ii<Tissue.Nb;ii++){
      file1 << Tissue.xbglobal(0,ii) << ", ";
      file1 << Tissue.xbglobal(1,ii) << endl;
    }
    file1.flush();
    file2.flush();

    printf("%f/%f\n",t,t_max);

    t = t+Tissue.dt;

  }   // for loop_num
  file1.close();
  return 0;
}
