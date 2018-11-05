//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <vector>
#include "cell.hpp"
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
float dt;      // time step
float len;     // Initial cell radius
int   NumLoop; // number of steps
int   Numcells;// number of cells

int main() {

  ReadParams(Numg,Nb,dims,cen,Src,rho,mu,dt,len,NumLoop,Numcells);

  mat stoch_xb = mat(2,Nb*Numcells,arma::fill::zeros);
  ofstream file1;
  tissue Tissue = tissue(Numg,dims,Nb,Src);

  for (int ii=0;ii<Numcells;ii++){
    Tissue.AddCell(len,3*ii*len,0);
  }
  Tissue.UpdateSources();

  file1.open ("output/boundarypositions.txt", ofstream::out);
  // Write initial data to file //
  for(int row = 0 ; row < Nb ; row++){
    file1 << Tissue.Cells[0].xb(0,row) << ", ";
    file1 << Tissue.Cells[0].xb(1,row) << endl;
    file1 << Tissue.Cells[1].xb(0,row) << ", ";
    file1 << Tissue.Cells[1].xb(1,row) << endl;
  }
  file1.flush();

  for (int loop_num=0; loop_num<NumLoop; loop_num++) {
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
    NavierStokes(Tissue.vg,Tissue.ug,Tissue.fg,Tissue.sg,Tissue.Ng,rho,mu,dt,Tissue.hg);
    Tissue.ug = Tissue.vg;

    //-- boundary velocities --//
    GridToBound(Tissue.ubglobal,Tissue.xbglobal,Tissue.Nb,Tissue.vg,Tissue.Ng,Tissue.hg,Tissue.hg,Tissue.xmin,Tissue.xmax);

    //-- new position of boundary points --//
    stoch_xb.randn();
    Tissue.xbglobal = Tissue.xbglobal + dt*Tissue.ubglobal + stoch_xb/1000;

    for (int ii=0;ii<Tissue.Nc;ii++){
      Tissue.Cells[ii].UpdateCom();
    }
    Tissue.Cells[1].UpdateCom();

    GlobalToLocal(Tissue);

    Tissue.UpdateSources();

    // Write data to file //
    for (int ii=0;ii<Tissue.Nc;ii++){
      for(int row = 0 ; row < Nb ; row++){
        file1 << Tissue.Cells[ii].xb(0,row) << ", ";
        file1 << Tissue.Cells[ii].xb(1,row) << endl;
      }
    }
    file1.flush();

    printf("%d/%d\n",loop_num+1,NumLoop);

  }   // for loop_num
  file1.close();
  return 0;
}
