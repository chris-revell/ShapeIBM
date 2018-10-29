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


using namespace std;
using namespace arma;

// System parameters
const int     Numg=16*16;                   // fluid grid size
const int     Nb=64;                        // number of boundary points
const int     dims=10;
const float   cen=0;                        // Fluid centre point
const float   Src=1;                      // source strength
const float   rho=1;                        // fluid density
const float   mu=1;                         // fluid viscosity
const float   dt=1;                      // time step
const float   len=0.2;                      // Initial cell radius
const int     NumLoop=40;                   // number of steps
const int     Nbs=2;                        // Number of fluid sources
const int     Nc=2;                         // Number of cells




int main() {

    mat sb   = mat(2,2,fill::zeros);
    mat ub   = mat(2,Nb,fill::zeros);     // Velocity of all boundary elements
//    mat fb   = mat(2,Nb,fill::zeros);     // Boundary forces
//    mat fadj = mat(2,Nb,fill::zeros);     // Adjacent forces array
//    mat fsec = mat(2,Nb,fill::zeros);     // Secondary forces array
//    mat fcen = mat(2,Nb,fill::zeros);     // Centre forces array
//    mat fopp = mat(2,Nb,fill::zeros);     // Opposite forces array


    cube fg = cube(Numg+1,Numg+1,2,fill::zeros);// Grid forces
    cube vg = cube(Numg+1,Numg+1,2,fill::zeros);// Grid velocities
    cube ug = cube(Numg+1,Numg+1,2,fill::zeros);// Previous grid velocities
    mat sbb= mat(2,2,fill::zeros);        //
    cube fvg0     = cube(Numg,Numg,2,fill::zeros);
    cube fvg1     = cube(Numg,Numg,2,fill::zeros);
    ofstream file1;
    tissue Tissue = tissue(Numg,dims);

    Tissue.AddCell(Nb,len,0,0);

    //-- distributed sources and sinks --//
    sb(0,0)  =cen;
    sb(1,0)  =cen;
    sb(0,1)  =-dims;
    sb(1,1)  =-dims;
    sbb(0,0) = Src;     // a source at the center
    sbb(0,1) =-Src;     // a sink in the corner

    file1.open ("output/boundarypositions.txt", std::ofstream::out | std::ofstream::app);
    // Write initial data to file //
    for(int row = 0 ; row < Nb ; row++){
        file1 << Tissue.Cells[0].xb(0,row) << ", ";
        file1 << Tissue.Cells[0].xb(1,row) << endl;
    }
    //file1 << "" << endl;
    file1.flush();

    for (int loop_num=0; loop_num<NumLoop; loop_num++) {
        //-- boundary forces --//
        Tissue.Cells[0].AdjacentForces();


        //-- grid sources --//
        BoundToGrid1(Tissue,sb,sbb,Nbs);
        //-- grid forces --//
        BoundToGrid2(fg,Tissue.Cells[0].xb,Tissue.Cells[0].fb,Tissue.Cells[0].Nb,Tissue.Ng,Tissue.hg,Tissue.hg,0.5*Tissue.hg,Tissue.xmin,Tissue.xmax);
        //-- compute grid velocity from NavierStokes --//
        NavierStokes(vg,ug,fg,Tissue.sg,Tissue.Ng,rho,mu,dt,Tissue.hg);
        ug = vg;
        //-- boundary velocities --//
        GridToBound(ub,Tissue.Cells[0].xb,Nb,vg,Tissue.Ng,Tissue.hg,Tissue.hg,Tissue.xmin,Tissue.xmax);
        //-- new position of boundary points --//
        ub.print();
        Tissue.Cells[0].xb = Tissue.Cells[0].xb + dt*ub;

        // Write data to file //
        for(int row = 0 ; row < Nb ; row++){
            file1 << Tissue.Cells[0].xb(0,row) << ", ";
            file1 << Tissue.Cells[0].xb(1,row) << endl;
        }

        //file1 << "" << endl;
        file1.flush();
        printf("%d/%d\n",loop_num+1,NumLoop);

    }   // for loop_num
    file1.close();
    return 0;
}
