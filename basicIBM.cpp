//C++ adaptation of CellGrowth.m from Katarzyna Rejniak
//#include <stdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

const int Ng=64;                  // fluid grid size
const int Nb=64;                  // number of boundary points
const float pi = 3.1415;

//--------------------------------------------------------------------//
// Define cell shape. Nb defines
// the number of boundary points; len defines a radius for a circular //
// cell; or a length of a side of the square                          //
//--------------------------------------------------------------------//
float DefineCellShape(const int& Nb, const float& len, const int& x, const int& y) {
  float hb=2*pi/static_cast<float>(Nb);
  if (x==1) {
    return (len*cos(y*hb));
  }
  else {
    return (len*sin(y*hb));
  }
}    // function DefineCellShape
//-----------------------------------------------------------------//

//-------------------------------------------------------------------//
// computes shifting of the index nn in the square periodic domain   //
// of size Ng+1 where the 1st and the (Ng+1)st columns and rows are  //
// identical according to the value of which:                        //
//    which = 1 -- nn shifted one element to the right               //
//    which =-1 -- nn shifted one element to the left                //
//-------------------------------------------------------------------//
int PeriodInd(const int& nn,const int& Ng,const int& which){
  int ind=0;
  if (which == (-1)){
    if (nn == 0){
      ind=Ng-1;
    }else{
      ind=nn-1;
    }
  }
  if (which == 1){
    if (nn == Ng){
      ind=1;
    }else{
      ind=nn+1;
    }
  }
  return ind;
} // function PeriodInd
//------------------------------------------------------------------//

//--------------------------------------------------------------------//
// transforms the real coordinate xy of the body into the correspon-  //
// ding coordinate pom inside the periodic domain (xmn,xmx)x(xmn,xmx) //                             //
//--------------------------------------------------------------------//
float IntoDom(const float& xy,const float& xmin,const float& xmax){
  float len;
  float pom;
  len=xmax-xmin;
  cout << "len "<< len << endl;
  pom=xy;
  cout << "pom "<< pom << endl;
  while (pom>xmax){
    pom=pom-len;
    cout << "while pom "<< pom << endl;
  }
  while (pom<xmin){
    pom=pom+len;
    cout << "while pom "<< pom << endl;
  }
  return pom;
} // function IntoDom
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// determines value of the unitary bell-shaped discrete approximation //
// to the Dirac delta function of radius h for a point located at    //
// distance r from the center.                                       //
//-------------------------------------------------------------------//
float DeltaFun(const float& r,const float& h){
  float dist;

  if (fabs(r) < (2*h)){
    dist=0.25*(1+cos(0.5*pi*r/h))/h;
  }else{
    dist=0;
  }
  return dist;
} // function DeltaFun
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// determines indices in1 and in2 of grid elements which need  to be //
// updated in the interpolation of point ll (shifted by ij according //
// to the reference point Nij). Depending on the location of point   //
// ll, there may be at most two active indeces for each element, due //
// to periodicity of fluid domain.                                   //
//   0 elements to update, if ll inside the domain,                  //
//   1 element  to update, if ll is close to the boundary,           //
//   2 elements to update, if ll is close to the corner.             //
//-------------------------------------------------------------------//
void IndDel(int& in1,int& in2,const float& ll,const int& ij,const int& Nij,const int& Nmx,const float& xmn,const float& xmx){
  // passive value - do nothing
  int pas=-100;

  in2=pas;
  if (ll < xmn){
    in1=Nmx+ij;
  }else if ((ll == xmn) || (ll == xmx)){
    in1=Nmx;
    in2=0;
  }else if (ll > xmx){
    in1=ij-1;
  }else{
    in1=Nij+ij-1;
  }
  return;
} // function IndDel
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//
// define adjacent forces between adjacent boundary points:           //
// (n-th to n-1 and n+1)                                              //
// Hookean springs between each adjacent boundary point with          //
// spring constant Spr and equilibrium distance Lrest=hb              //
//--------------------------------------------------------------------//
void AdjacentForces(mat& fbb,const mat& xb,const int& Nb,const float& hb, const float& Spr) {
  float Lrest;
  float dl1;
  float dl2;
  float dr1;
  float dr2;
  float ndl;
  float ndr;
  for (int ii=0; ii<Nb; ii++) {
    Lrest=hb;
    if (ii==0) {
      dl1=xb(0,Nb-1)-xb(0,ii);
      dl2=xb(1,Nb-1)-xb(1,ii);
      dr1=xb(0,ii+1)-xb(0,ii);
      dr2=xb(1,ii+1)-xb(1,ii);
    } else if (ii==Nb-1) {
      dl1=xb(0,ii-1)-xb(0,ii);
      dl2=xb(1,ii-1)-xb(1,ii);
      dr1=xb(0,0)-xb(0,ii);
      dr2=xb(1,0)-xb(1,ii);
    } else {
      dl1=xb(0,ii-1)-xb(0,ii);
      dl2=xb(1,ii-1)-xb(1,ii);
      dr1=xb(0,ii+1)-xb(0,ii);
      dr2=xb(1,ii+1)-xb(1,ii);
    }
    ndl=sqrt(pow(dl1,2)+pow(dl2,2));
    ndr=sqrt(pow(dr1,2)+pow(dr2,2));
    fbb(0,ii)=Spr*(ndl-Lrest)*dl1/ndl+Spr*(ndr-Lrest)*dr1/ndr;
    fbb(1,ii)=Spr*(ndl-Lrest)*dl2/ndl+Spr*(ndr-Lrest)*dr2/ndr;
  }
  return;
}// function AdjacentForces
//----------------------------------------------------------------------//

//----------------------------------------------------------------------//
// define secondary adjacent forces between second to adjacent boundary //
// points: (n-th to n-2 and n+2)                                        //
// Hookean springs with spring constant Spr again but this time         //
// equilibrium distance Lrest=2*hb                                      //
//----------------------------------------------------------------------//
void SecondaryForces(mat& fbb,const mat& xb,const int& Nb,const float& hb,const float& Spr,const int& connect) {
  float dl1;
  float dl2;
  float dr1;
  float dr2;
  float ndl;
  float ndr;
  float Lrest;
  if (connect==1) {
    Lrest=2*hb;
    for (int ii=0; ii<Nb; ii++) {
      if (ii==0) {
        dl1=xb(0,Nb-2)-xb(0,ii);
        dl2=xb(1,Nb-2)-xb(1,ii);
        dr1=xb(0,ii+2)-xb(0,ii);
        dr2=xb(1,ii+2)-xb(1,ii);
      } else if (ii==1) {
        dl1=xb(0,Nb-1)-xb(0,ii);
        dl2=xb(1,Nb-1)-xb(1,ii);
        dr1=xb(0,ii+2)-xb(0,ii);
        dr2=xb(1,ii+2)-xb(1,ii);
      } else if (ii==Nb-2) {
        dl1=xb(0,ii-2)-xb(0,ii);
        dl2=xb(1,ii-2)-xb(1,ii);
        dr1=xb(0,0)-xb(0,ii);
        dr2=xb(1,0)-xb(1,ii);
      } else if (ii==Nb-1) {
        dl1=xb(0,ii-2)-xb(0,ii);
        dl2=xb(1,ii-2)-xb(1,ii);
        dr1=xb(0,1)-xb(0,ii);
        dr2=xb(1,1)-xb(1,ii);
      } else {
        dl1=xb(0,ii-2)-xb(0,ii);
        dl2=xb(1,ii-2)-xb(1,ii);
        dr1=xb(0,ii+2)-xb(0,ii);
        dr2=xb(1,ii+2)-xb(1,ii);
      }
      ndl=sqrt(pow(dl1,2)+pow(dl2,2));
      ndr=sqrt(pow(dr1,2)+pow(dr2,2));
      fbb(0,ii)=Spr*(ndl-Lrest)*dl1/ndl+Spr*(ndr-Lrest)*dr1/ndr;
      fbb(1,ii)=Spr*(ndl-Lrest)*dl2/ndl+Spr*(ndr-Lrest)*dr2/ndr;
    }
  }
  return;
} // SecondaryForces
//--------------------------------------------------------------------//

//----------------------------------------------------------------------//
// define central forces between boundary points and the cell nucleus   //
// only half of all boundary points are connected to cell nucleus to    //
// allow the whole cell to grow                                         //
// Again Hookean forces with spring constant Spr                        //
// and equilibium radius len where len is the typical cell radius       //
//----------------------------------------------------------------------//
void CenterForces(mat& fbb,const mat& xb,const int& Nb,const float& cen,const float& len,const float& Spr,const int& connect){
  float dl1;
  float dl2;
  float ndl;
  float Lrest;

  if (connect==2){
    Lrest=len;
    for (int ii=0; ii<Nb/2; ii++){
      dl1=cen-xb(0,ii);
      dl2=cen-xb(1,ii);
      ndl=sqrt(pow(dl1,2)+pow(dl2,2));

      fbb(0,ii)=Spr*(ndl-Lrest)*dl1/ndl;
      fbb(1,ii)=Spr*(ndl-Lrest)*dl2/ndl;
    }
  }
  return;
} // function CenterForces
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//
// define opposite forces between opposite boundary points            //
// Again Hookean springs with spring constant Spr and this time       //
// equilibrium radius Lrest=2*len where len is the typical cell radius//
//--------------------------------------------------------------------//
void OppositeForces(mat& fbb,const mat& xb,const int& Nb,const float& len,const float& Spr,const int& connect){
  float Lrest;
  float dl1;
  float dl2;
  float dr1;
  float dr2;
  float ndl;
  float ndr;
  int Nb2;

  if (connect==3){
    Lrest=2*len;
    Nb2=floor(Nb/4);
    for (int ii=0; ii<Nb2; ii++){
      dl1=xb(0,3*Nb2-ii)-xb(0,ii);
      dl2=xb(1,3*Nb2-ii)-xb(1,ii);
      ndl = sqrt(pow(dl1,2)+pow(dl2,2));
      fbb(0,ii)       = Spr*(ndl-Lrest)*dl1/ndl;
      fbb(1,ii)       = Spr*(ndl-Lrest)*dl2/ndl;
      fbb(0,3*Nb2-ii) = -Spr*(ndl-Lrest)*dl1/ndl;
      fbb(1,3*Nb2-ii) = -Spr*(ndl-Lrest)*dl2/ndl;
    }
  }
  return;
}// function OppositeForces
//--------------------------------------------------------------------//

//-----------------------------------------------------------------//
// spreads the material values sb(Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1) in the //
// square domain (xmin,xmax)^2 with mesh width hg, material points //
// separation hb and a radius of the discrete delta function hdl.  //
//-----------------------------------------------------------------//
void BoundToGrid1(mat& sg,const mat& xb,const mat& sb,const int& Nb,const int& Ng,const float& hdl,const float& hg,const float& hb,const float& xmn,const float& xmx){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  float xbb0;
  float xbb1;
  int Nx,Ny;

  for (int n3=0;n3<Nb;n3++){

    // Move points into the domain. xbb0 and xbb1 are the coordinates of the point in the square domain.
    xbb0=IntoDom(xb(0,n3),xmn,xmx);
    xbb1=IntoDom(xb(1,n3),xmn,xmx);
    //cout << "xbb0 " << xbb0 << endl;
    //cout << "xbb1 " << xbb1 << endl;
    // determine indices of the nearest lower-down grid point
    Nx=1+floor((xbb0-xmn)/hg);
    Ny=1+floor((xbb1-xmn)/hg);
    //cout << "Nx " << Nx << endl;
    //cout << "Ny " << Ny << endl;
    // tests all 16 possible grid points
    for (int ii=-1; ii<3; ii++){
      for (int jj=-1; jj<3; jj++){
        //cout << "ii " << ii << endl;
        //cout << "jj " << jj << endl;
        // compute the interpolation Delta function
        llx=xmn+(Nx-1)*hg+ii*hg;
        rr=fabs(xbb0-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmn+(Ny-1)*hg+jj*hg;
        rr=fabs(xbb1-lly);
        dy=DeltaFun(rr,hdl);
        //cout << "dx " << dx << endl;
        //cout << "dy " << dy << endl;
        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmn,xmx);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmn,xmx);
        //cout << "x1 " << x1 << endl;
        //cout << "x2 " << x2 << endl;
        //cout << "y1 " << y1 << endl;
        //cout << "y2 " << y2 << endl;
        // update the values if points are not passive
        if (dx*dy > 0){
          sg(x1,y1)  = sg(x1,y1) + sb(0,n3)*dx*dy*hb;
          if (x2 != pas){
            sg(x2,y1)= sg(x2,y1) + sb(0,n3)*dx*dy*hb;
          }
          if (y2 != pas){
            sg(x1,y2)= sg(x1,y2) + sb(0,n3)*dx*dy*hb;
          }
          if ((x2 != pas) & (y2 != pas)){
            sg(x2,y2)= sg(x2,y2) + sb(0,n3)*dx*dy*hb;
          }
        }
      }  // for jj
    } // for ii
  } // for n3
  return;
} // function BoundToGrid1
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// spreads the material values sb(1,Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1,2) in the //
// square domain (xmin,xmax)^2 with mesh width hg, material points   //
// separation hb and a radius of the discrete delta function hdl.    //
//-------------------------------------------------------------------//
void BoundToGrid2(cube& sg,const mat& xb,const mat& sb,const int& Nb,const int& Ng,const float& hdl,const float& hg,const float& hb,const float& xmin,const float& xmax){

  // passive value - do nothing
  int pas=-100;
  int Nx,Ny;
  float xbb0,xbb1;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  float llx,rr,dx,lly,dy;

  for (int n3=0; n3<Nb; n3++){
    // move points into the domain

    xbb0=IntoDom(xb(0,n3),xmin,xmax);

    xbb1=IntoDom(xb(1,n3),xmin,xmax);

    // determine indices of the nearest lower-down grid point
    Nx=1+floor((xbb0-xmin)/hg);
    Ny=1+floor((xbb1-xmin)/hg);

    // tests all 16 possible grid points
    for (int ii=-1; ii<2;ii++){
      for (int jj=-1;jj<2;jj++){
        // compute the interpolation Delta function
        llx=xmin+(Nx-1)*hg+ii*hg;
        rr= fabs(xbb0-llx);
        dx= DeltaFun(rr,hdl);
        lly=xmin+(Ny-1)*hg+jj*hg;
        rr= fabs(xbb1-lly);
        dy= DeltaFun(rr,hdl);

        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);

        // update the values if points are not pasive
        if (dx*dy > 0){
          sg(x1,y1,0)=sg(x1,y1,0)+sb(0,n3)*dx*dy*hb;
          sg(x1,y1,1)=sg(x1,y1,1)+sb(1,n3)*dx*dy*hb;

          if (x2 != pas){
            sg(x2,y1,0)=sg(x2,y1,0)+sb(0,n3)*dx*dy*hb;
            sg(x2,y1,1)=sg(x2,y1,1)+sb(1,n3)*dx*dy*hb;
          }
          if (y2 != pas){
            sg(x1,y2,0)=sg(x1,y2,0)+sb(0,n3)*dx*dy*hb;
            sg(x1,y2,1)=sg(x1,y2,1)+sb(1,n3)*dx*dy*hb;
          }
          if ((x2 != pas) && (y2 != pas)){
            sg(x2,y2,0)=sg(x2,y2,0)+sb(0,n3)*dx*dy*hb;
            sg(x2,y2,1)=sg(x2,y2,1)+sb(1,n3)*dx*dy*hb;
          }
        }

      } // jj
    } // ii
  } // for n3
  return;
} // function BoundToGrid2
//---------------------------------------------------------------------//

//---------------------------------------------------------------------//
// Uses the Fast Fourier Method to solve the Navier-Stokes equations   //
// for updating grid velocity ug due to grid forces fg and grid source //
// distribution sg, rho and mu are fluid constants, dt is a time step, //
// hg is a mesh width.                                                 //
//---------------------------------------------------------------------//
void NavierStokes(cube& vg,const cube& ug,const cube& fg,const mat& sg,const int& Ng,const float& rho,const float& mu,const float& dt,const float& hg){

  float pom,B1,B2,Aa,Bb,Bv;
  int in1,in2;
  float Eps=0.0000001;
  mat dummymat;
  cx_mat fvgg     = cx_mat(Ng,Ng,fill::zeros);
  cx_mat fsg      = cx_mat(Ng,Ng,fill::zeros);
  cx_mat fug1     = cx_mat(Ng,Ng,fill::zeros);
  cx_mat fug2     = cx_mat(Ng,Ng,fill::zeros);
  cube fvg0     = cube(Ng,Ng,2,fill::zeros);
  cube fvg1     = cube(Ng,Ng,2,fill::zeros);
  cx_mat vg0      = cx_mat(Ng,Ng,fill::zeros);
  cx_mat vg1      = cx_mat(Ng,Ng,fill::zeros);

  // stage n terms: force density fg, source distribution sg and current
  // velocity ug
  for (int n1=0; n1<Ng+1; n1++){
    for (int n2=0; n2<Ng+1; n2++){
      for (int ik=0; ik<2; ik++){
        // upwind scheme for the advection term
        if (ug(n1,n2,0) < 0){
          in1=PeriodInd(n1,Ng,1);
          //cout << "n1,n2,ik " << endl;
          //cout << n1 << endl;
          //cout << n2 << endl;
          //cout << ik << endl;
          //cout << "in1 " << in1 << endl;
          pom=ug(in1,n2,ik)-ug(n1,n2,ik);
          //cout << "pom " << pom << endl;
        }else{
          in1=PeriodInd(n1,Ng,-1);
          //cout << "n1,n2,ik " << endl;
          //cout << n1 << endl;
          //cout << n2 << endl;
          //cout << ik << endl;
          //cout << "in1 " << in1 << endl;
          pom=ug(n1,n2,ik)-ug(in1,n2,ik);
          //cout << "pom " << pom << endl;
        }
        vg(n1,n2,ik)=ug(n1,n2,0)*pom;

        if (ug(n1,n2,1) < 0){
          in2=PeriodInd(n2,Ng,1);
          pom=ug(n1,in2,ik)-ug(n1,n2,ik);
        }else{
          in2=PeriodInd(n2,Ng,-1);
          pom=ug(n1,n2,ik)-ug(n1,in2,ik);
        }
        //cout << "n1,n2,ik " << endl;
        //cout << n1 << endl;
        //cout << n2 << endl;
        //cout << ik << endl;
        //cout << "vg(n1,n2,ik) " << vg(n1,n2,ik) << endl;
        //cout << "pom " << pom << endl;
        vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,1)*pom;
        vg(n1,n2,ik)=-dt*vg(n1,n2,ik)/hg;
        //cout << "vg(n1,n2,ik) " << vg(n1,n2,ik) << endl;

        // central difference for the grad of source term
        if (ik == 0){
          in1=PeriodInd(n1,Ng,1);
          in2=PeriodInd(n1,Ng,-1);
          pom=sg(in1,n2)-sg(in2,n2);
        }else if (ik == 1){
          in1=PeriodInd(n2,Ng,1);
          in2=PeriodInd(n2,Ng,-1);
          pom=sg(n1,in1)-sg(n1,in2);
        }
        vg(n1,n2,ik)=vg(n1,n2,ik)+dt*mu*pom/(6*hg*rho*rho);
        //cout << "vg(n1,n2,ik) " << vg(n1,n2,ik) << endl;
        // current vlocity and force terms
        vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,ik)+dt*fg(n1,n2,ik)/rho;
        //cout << "vg(n1,n2,ik) " << vg(n1,n2,ik) << endl;
      } // for ik
    } // for n2
  } // for n1

  // the Fast Fourier transforms of source distribution sg and stage n term vg
  fsg  = fft(sg(span(0,Ng-1),span(0,Ng-1)));
  //cout << sg(span(0,Ng-1),span(0,Ng-1))(0,0) << endl;
  //cout << "fsg(0,0) " << fsg(0,0) << endl;
  dummymat = vg.slice(0);
  //cout << dummymat(0,0) << endl;
  fug1 = fft(dummymat(span(0,Ng-1),span(0,Ng-1)));
  //cout << "fug1(0,0) " << fug1(0,0) << endl;
  dummymat = vg.slice(1);
  //cout << dummymat(0,0) << endl;
  fug2 = fft(dummymat(span(0,Ng-1),span(0,Ng-1)));
  //cout << "fug2(0,0) " << fug2(0,0) << endl;
  // determines fug - the Fourier Transform of the velocity field at the stage n+1
  for (int n1=0;n1<Ng;n1++){
    for (int n2=0;n2<Ng;n2++){
      B1=sin(2*pi*n1/Ng);
      B2=sin(2*pi*n2/Ng);
      Bb=pow(B1,2)+pow(B2,2);
      Aa=1+4*mu*dt*(pow(sin(pi*n1/Ng),2)+pow(sin(pi*n2/Ng),2))/(rho*hg*hg);

      B1=sin(2*pi*n1/Ng);
      B2=sin(2*pi*n2/Ng);
      Bb=pow(B1,2)+pow(B2,2);
      Aa=1+4*mu*dt*(pow(sin(pi*n1/Ng),2)+pow(sin(pi*n2/Ng),2))/(rho*hg*hg);

      if (Bb < Eps){
        fvg0(n1,n2,0)=real(fug1(n1,n2))/Aa;
        fvg0(n1,n2,1)=imag(fug1(n1,n2))/Aa;
        fvg1(n1,n2,0)=real(fug2(n1,n2))/Aa;
        fvg1(n1,n2,1)=imag(fug2(n1,n2))/Aa;
      }else{
        Bv=B1*real(fug1(n1,n2))+B2*real(fug2(n1,n2));
        fvg0(n1,n2,0)=(Bb*real(fug1(n1,n2))-B1*Bv)/(Aa*Bb);
        fvg1(n1,n2,0)=(Bb*real(fug2(n1,n2))-B2*Bv)/(Aa*Bb);

        Bv=B1*imag(fug1(n1,n2))+B2*imag(fug2(n1,n2));
        fvg0(n1,n2,1)=(Bb*imag(fug1(n1,n2))-B1*Bv)/(Aa*Bb);
        fvg1(n1,n2,1)=(Bb*imag(fug2(n1,n2))-B2*Bv)/(Aa*Bb);

        fvg0(n1,n2,0)=fvg0(n1,n2,0)+hg*B1*imag(fsg(n1,n2))/(Bb*rho);
        fvg0(n1,n2,1)=fvg0(n1,n2,1)-hg*B1*real(fsg(n1,n2))/(Bb*rho);
        fvg1(n1,n2,0)=fvg1(n1,n2,0)+hg*B2*imag(fsg(n1,n2))/(Bb*rho);
        fvg1(n1,n2,1)=fvg1(n1,n2,1)-hg*B2*real(fsg(n1,n2))/(Bb*rho);
      }
    } // for n2
  } // for n1

  // the inverse Fast Fourier Method of fvg
  dummymat = fvg0.slice(0);
  fvgg.set_real(dummymat(span(0,Ng-1),span(0,Ng-1)));
  dummymat = fvg0.slice(1);
  fvgg.set_imag(dummymat(span(0,Ng-1),span(0,Ng-1)));
  vg0 = ifft(fvgg);

  dummymat = fvg1.slice(0);
  fvgg.set_real(dummymat(span(0,Ng-1),span(0,Ng-1)));
  dummymat = fvg1.slice(1);
  fvgg.set_imag(dummymat(span(0,Ng-1),span(0,Ng-1)));
  vg1 = ifft(fvgg);

  for (int ii=0; ii<Ng; ii++){
    for (int jj=0; jj<Ng; jj++){
      vg(ii,jj,0)=real(vg0(ii,jj));
      vg(ii,jj,1)=real(vg1(ii,jj));
    }
  }

  for (int ii=0; ii<Ng; ii++){
    vg(Ng,ii,0)=vg(0,ii,0);
    vg(ii,Ng,0)=vg(ii,Ng,0);
    vg(Ng,ii,1)=vg(0,ii,1);
    vg(ii,Ng,1)=vg(ii,Ng,1);
  }
} // function NavierStokes
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
// interpolates the grid values fg(Ng,Ng,2) (velocities) to the material //
// values fb(1,Nb) defined at the material points xb(1,Nb) in the square //
// domain (xmn,xmx)^2 with mesh width hg and a radius of the discrete    //
// Dirac delta hdl.                                                      //
//-----------------------------------------------------------------------//
void GridToBound(mat& fb,const mat& xb,const int& Nb,const cube& fg,const int& Ng,const float& hdl,const float& hg,const float& xmn,const float& xmx){

  float llx,rr,dx,lly,dy,xbb0,xbb1;
  int Nx,Ny;
  int x1 = 0;
  int x2 = 0;
  int y1 = 0;
  int y2 = 0;


  for (int n3=0; n3<1; n3++){
    // moves points into the domain
    cout << "" << endl;
    cout << "" << endl;
    cout << "" << endl;
    //cout << "xb(0,n3)" << xb(0,n3) << endl;
    //cout << "xb(1,n3)" << xb(1,n3) << endl;
    xbb0=IntoDom(xb(0,n3),xmn,xmx);
    xbb1=IntoDom(xb(1,n3),xmn,xmx);
    //cout << "xbb0" << xbb0 << endl;
    //cout << "xbb1" << xbb1 << endl;
    //computes the indices of the nearest down-left grid point
    Nx=1+floor((xbb0-xmn)/hg);
    Ny=1+floor((xbb1-xmn)/hg);
    //cout << "Nx" << Nx << endl;
    //cout << "Ny" << Ny << endl;
    // test all 16 neighboring points
    for (int ii=-1;ii<3;ii++){
      for (int jj=-1;jj<3;jj++){
        // determine the value of the interpolation Delta-function
        //cout << "" << endl;
        //cout << "ii,jj" << ii << jj << endl;
        llx=xmn+(Nx-1)*hg+ii*hg;
        rr=fabs(xbb0-llx);
        //cout << "llx" << llx << endl;
        //cout << "rr" << rr << endl;
        dx=DeltaFun(rr,hdl);
        lly=xmn+(Ny-1)*hg+jj*hg;
        rr=fabs(xbb1-lly);
        dy=DeltaFun(rr,hdl);
        //cout << "dx" << dx << endl;
        //cout << "dy" << dy << endl;
        // determine the indices of grid points gaining positive impact
        IndDel(x1,x2,llx,ii,Nx,Ng,xmn,xmx);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmn,xmx);
        //cout << "x1" << x1 << endl;
        //cout << "x2" << x2 << endl;
        //cout << "y1" << y1 << endl;
        //cout << "y2" << y2 << endl;
        // update the values if inside the impact domain
        if (dx*dy > 0){
          fb(0,n3)=fb(0,n3)+fg(x1,y1,0)*dx*dy*hg*hg;
          fb(1,n3)=fb(1,n3)+fg(x1,y1,1)*dx*dy*hg*hg;
        }
        cout << "fg(x1,y1,0)" << fg(x1,y1,0) << endl;
        cout << "fg(x1,y1,1)" << fg(x1,y1,1) << endl;
        cout << "dx" << dx << endl;
        cout << "dy" << dy << endl;
        cout << "hg" << hg << endl;
        cout << "fb(0,n3)" << fb(0,n3) << endl;
        cout << "fb(1,n3)" << fb(1,n3) << endl;

      } // for jj
    } // for ii
  }// for n3
  return;
} // function GridToBound
//-------------------------------------------------------------------//

int main() {

  float xmin=-1;
  float xmax=1;
  float cen=(xmax+xmin)/2;      // fluid domain (square)
  float hg=float(xmax-xmin)/float(Ng);    // fluid mesh width
  float Src=1;                // source strength
  float Spr=100;              // spring stiffness
  float rho=1;                // fluid density
  float mu=1;                 // fluid viscosity
  float dt=0.05;              // time step
  float len=0.2;              // body radius/half length
  int connect=3;              // spring connections:
                              // 0-adjacent only, 1-adjacent and secondary
                              // 2-adjacent and center, 3-adjcent and opposite
  float scaleF=1.0;           // scaling parameter for forces
  float scaleV=1.0;           // scaling parameter for velocities
  int NumLoop=1;             // number of steps
  int mod_num=5;              // frequency
  int Nbs=2;
  float hb;                   // Spacing between boundary elements?
  mat sb   = mat(2,2,fill::zeros);
  mat xb   = mat(2,Nb,fill::zeros);     // Position of all boundary elements
  mat ub   = mat(2,Nb,fill::zeros);     // Velocity of all boundary elements
  mat fb   = mat(2,Nb,fill::zeros);     // Boundary forces
  mat fadj = mat(2,Nb,fill::zeros);     // Adjacent forces array
  mat fsec = mat(2,Nb,fill::zeros);     // Secondary forces array
  mat fcen = mat(2,Nb,fill::zeros);     // Centre forces array
  mat fopp = mat(2,Nb,fill::zeros);     // Opposite forces array
  cube xg = cube(Ng+1,Ng+1,2,fill::zeros);// Fluid grid array
  mat sg = mat(Ng+1,Ng+1,fill::zeros);  // Fluid grid
  cube fg = cube(Ng+1,Ng+1,2,fill::zeros);// Grid forces
  cube vg = cube(Ng+1,Ng+1,2,fill::zeros);// Grid velocities
  cube ug = cube(Ng+1,Ng+1,2,fill::zeros);// Previous grid velocities
  mat sbb= mat(2,2,fill::zeros);        //

  cube fvg0     = cube(Ng,Ng,2,fill::zeros);
  cube fvg1     = cube(Ng,Ng,2,fill::zeros);

  //-- define fluid grid --//
  for (int ii=0;ii<Ng+1;ii++){
    for (int jj=0;jj<Ng+1;jj++){
      xg(ii,jj,0)=xmin+ii*hg;
      xg(ii,jj,1)=xmin+jj*hg;
    }
  }

  //-- cell shape --//
  for (int ii=0;ii<Nb;ii++) {
    xb(0,ii) = DefineCellShape(Nb,len,1,ii);
    xb(1,ii) = DefineCellShape(Nb,len,2,ii);
  }

  hb=sqrt(pow((xb(0,1)-xb(0,2)),2)+pow((xb(1,1)-xb(1,2)),2));

  //-- distributed sources and sinks --//
  sb(0,0)  =cen;
  sb(1,0)  =cen;
  sb(0,1)  =xmin;
  sb(1,1)  =xmin;
  sbb(0,0) = Src;     // a source at the center
  sbb(0,1) =-Src;     // a sink in the corner

  ofstream file1 ("boundarypositions.txt");

  for (int loop_num=0; loop_num<NumLoop; loop_num++) {
    //-- boundary forces --//
    //cout << "xb(0,0) " << xb(0,0) << endl;
    //cout << "xb(0,Nb-1) " << xb(0,Nb-1) << endl;
    AdjacentForces(fadj,xb,Nb,hb,Spr);            // adjacent
    //cout << "fadj(0,0) " << fadj(0,0) << endl;
    SecondaryForces(fsec,xb,Nb,hb,Spr,connect);   // secondary
    //cout << "fsec(0,0) " << fsec(0,0) << endl;
    CenterForces(fcen,xb,Nb,cen,len,Spr,connect); // center
    //cout << "fcen(0,0) " << fcen(0,0) << endl;
    OppositeForces(fopp,xb,Nb,len,Spr,connect);   // opposite
    //cout << "fopp(0,0) " << fopp(0,0) << endl;
    for (int ii=0; ii<Nb; ii++){
      fb(0,ii)=fadj(0,ii)+fsec(0,ii)+fcen(0,ii)+fopp(0,ii);    // add all forces
      fb(1,ii)=fadj(1,ii)+fsec(1,ii)+fcen(1,ii)+fopp(1,ii);    // add all forces
    }
    //cout << "fb(0,0) " << fb(0,0) << endl;
    //-- grid sources --//
    BoundToGrid1(sg,sb,sbb,Nbs,Ng,hg,hg,0.5*hg,xmin,xmax);
    //cout << "sg(0,0) " << sg(0,0) << endl;
    //cout << "sg(0,1) " << sg(0,1) << endl;
    //cout << "sg(1,0) " << sg(1,0) << endl;
    //cout << "sg(1,1) " << sg(1,1) << endl;
    //-- grid forces --//
    BoundToGrid2(fg,xb,fb,Nb,Ng,hg,hg,0.5*hg,xmin,xmax);
    //cout << "fg(0,0,0) " << fg(0,0,0) << endl;
    //-- compute grid velocity from NavierStokes --//
    //cout << "ug(0,0,0) " << ug(0,0,0) << endl;
    NavierStokes(vg,ug,fg,sg,Ng,rho,mu,dt,hg);
    ug = vg;
    cout << "vg(0,0,0) " << vg(0,0,0) << endl;
    cout << "vg(0,0,1) " << vg(0,0,1) << endl;
    //-- boundary velocities --//
    GridToBound(ub,xb,Nb,vg,Ng,hg,hg,xmin,xmax);
    cout << "ub(0,0) " << ub(0,0) << endl;
    //-- new position of boundary points --//
    xb = xb + dt*fb;

    // Write data to file //
    for(int row = 0 ; row < Nb ; row++){
      file1 << xb(0,row) << ", ";
      file1 << xb(1,row) << endl;
    }
    file1 << "" << endl;

  }   // for loop_num
  file1.close();
  return 0;
}
