//C++ adaptation of CellGrowth.m from Katarzyna Rejniak

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
//#include<cmath>
//#include "mkl_dfti.h"
//#include "array"
//using std::array;
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
using namespace std;


const int Ng=64;                  // fluid grid size
const int Nb=64;                  // number of boundary points
const float pi = 3.1415;

void fourn(float data[], unsigned long nn[], int ndim, int isign) { //Replaces data by its ndim-dimensional discrete Fourier transform, if isign is input as 1. nn[1..ndim] is an integer array containing the lengths of each dimension (number of complex values), which MUST all be powers of 2. data is a real array of length twice the product of these lengths, in which the data are stored as in a multidimensional complex array: real and imaginary parts of each element are in consecutive locations, and the rightmost index of the array increases most rapidly as one proceeds along data. For a two-dimensional array, this is equivalent to storing the array by rows. If isign is input as −1, data is replaced by its inverse transform times the product of the lengths of all dimensions.

  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;
  for (ntot=1,idim=1;idim<=ndim;idim++) ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim]; nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
        for (i1=i2;i1<=i2+ip1-2;i1+=2) {
          for (i3=i1;i3<=ip3;i3+=ip2) {
            i3rev=i2rev+i3-i2;
            SWAP(data[i3],data[i3rev]);
            SWAP(data[i3+1],data[i3rev+1]);
          }
        }
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
        for (i1=i3;i1<=i3+ip1-2;i1+=2){
          for (i2=i1;i2<=ip3;i2+=ifp2) {
            k1=i2;
            k2=k1+ifp1;
            tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
            tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
            data[k2]=data[k1]-tempr;
            data[k2+1]=data[k1+1]-tempi;
            data[k1] += tempr;
            data[k1+1] += tempi;
          }
        }
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
  nprev *= n;
  }
}

void rlft3(float ***data, float **speq, unsigned long nn1, unsigned long nn2, unsigned long nn3, int isign)
//Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3] (where nn1 = 1 for the case of a logically
//two-dimensional array), this routine returns (for isign=1) the complex fast Fourier transform as two complex arrays:
//On output, data contains the zero and positive frequency values of the third frequency component,
//while speq[1..nn1][1..2*nn2] contains the Nyquist critical frequency values of the third frequency component.
//First (and second) frequency components are stored for zero, positive, and negative frequencies,
//in standard wrap- around order. See text for description of how complex values are arranged. For isign=-1,
//the inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative factor) is performed,
//with output data (viewed as a real array) deriving from input data (viewed as complex) and speq.
//For inverse transforms on data not generated first by a forward transform, make sure the complex
//input data array satisfies property (12.5.2). The dimensions nn1, nn2, nn3 must always be integer powers of 2.
{
  void fourn(float data[], unsigned long nn[], int ndim, int isign);
  void nrerror(char error_text[]);
  unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
  double theta,wi,wpi,wpr,wr,wtemp;
  float c1,c2,h1r,h1i,h2r,h2i;
  //if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3) nrerror("rlft3: problem with dimensions or contiguity of data array\n");
  c1=0.5;
  c2 = -0.5*isign; theta=isign*(6.28318530717959/nn3); wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp; wpi=sin(theta);
  nn[1]=nn1;
  nn[2]=nn2;
  nn[3]=nn3 >> 1;
  if (isign == 1) {
    fourn(&data[1][1][1]-1,nn,3,isign);
    for (i1=1;i1<=nn1;i1++) {
      for (i2=1,j2=0;i2<=nn2;i2++) {
        speq[i1][++j2]=data[i1][i2][1];
        speq[i1][++j2]=data[i1][i2][2];
      }
    }
  }
  for (i1=1;i1<=nn1;i1++) {
    j1=(i1 != 1 ? nn1-i1+2 : 1);
    wr=1.0;
    wi=0.0;
    for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
      for (i2=1;i2<=nn2;i2++) {
        if (i3 == 1) {
          j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1); h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
          h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
          h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
          h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
          data[i1][i2][1]=h1r+h2r;
          data[i1][i2][2]=h1i+h2i;
          speq[j1][j2]=h1r-h2r;
          speq[j1][j2+1]=h2i-h1i;
        } else {
          j2=(i2 != 1 ? nn2-i2+2 : 1);
          j3=nn3+3-(i3<<1);
          h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
          h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
          h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
          h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
          data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
          data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
          data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
          data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
        }
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  }
  if (isign == -1)
  fourn(&data[1][1][1]-1,nn,3,isign);
}

//--------------------------------------------------------------------//
// Define cell shape. Nb defines
// the number of boundary points; len defines a radius for a circular //
// cell; or a length of a side of the square                          //
//--------------------------------------------------------------------//
int DefineCellShape(int Nb, float len, int x, int y) {
  float hb=2*pi/Nb;
  if (x==1) {
    return (len*cos((y-1)*hb));
  }
  else {
    return (len*sin((y-1)*hb));
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
int PeriodInd(int nn,int Ng,int which){
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
// ding coordinate pom inside the periodic domain [xmn,xmx]x[xmn,xmx] //                             //
//--------------------------------------------------------------------//
int IntoDom(float xy,int xmin,int xmax){
  int len;
  int pom;
  len=xmax-xmin;
  pom=xy;
  while (pom>xmax){
    pom=pom-len;
  }
  while (pom<xmin){
    pom=pom+len;
  }
  return pom;
} // function IntoDom
//-------------------------------------------------------------------//

//-------------------------------------------------------------------//
// determines value of the unitary bell-shaped discrete approximation //
// to the Dirac delta function of radius h for a point located at    //
// distance r from the center.                                       //
//-------------------------------------------------------------------//
float DeltaFun(float r,float h){
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
void IndDel(int in1,int in2,float ll,int ij,int Nij,int Nmx,int xmn,int xmx){
  // passive value - do nothing
  int pas=-100;

  in2=pas;
  if (ll < xmn){
    in1=Nmx+1+ij;
  }else if ((ll == xmn) || (ll == xmx)){
    in1=Nmx+1;
    in2=1;
  }else if (ll > xmx){
    in1=ij;
  }else{
    in1=Nij+ij;
  }
  return;
} // function IndDel
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//
// define adjacent forces between adjacent boundary points:           //
// (n-th to n-1 and n+1)                                              //
//--------------------------------------------------------------------//
void AdjacentForces(float fbb[2][Nb],float xb[2][Nb],int Nb,float hb, float Spr) {
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
      dl1=xb[0][Nb-1]-xb[0][ii];
      dl2=xb[1][Nb-1]-xb[1][ii];
      dr1=xb[0][ii+1]-xb[0][ii];
      dr2=xb[1][ii+1]-xb[1][ii];
    } else if (ii==Nb-1) {
      dl1=xb[0][ii-1]-xb[0][ii];
      dl2=xb[1][ii-1]-xb[1][ii];
      dr1=xb[0][0]-xb[0][ii];
      dr2=xb[1][0]-xb[1][ii];
    } else {
      dl1=xb[0][ii-1]-xb[0][ii];
      dl2=xb[1][ii-1]-xb[1][ii];
      dr1=xb[0][ii+1]-xb[0][ii];
      dr2=xb[1][ii+1]-xb[1][ii];
    }
    ndl=sqrt(pow(dl1,2)+pow(dl2,2));
    ndr=sqrt(pow(dr1,2)+pow(dr2,2));
    fbb[0][ii]=fbb[0][ii]+(Spr*(ndl-Lrest)*dl1/ndl)+(Spr*(ndr-Lrest)*dr1/ndr);
    fbb[1][ii]=fbb[1][ii]+(Spr*(ndl-Lrest)*dl2/ndl)+(Spr*(ndr-Lrest)*dr2/ndr);
  }
  return;
}// function AdjacentForces
//----------------------------------------------------------------------//

//----------------------------------------------------------------------//
// define secondary adjacent forces between second to adjacent boundary //
// points: (n-th to n-2 and n+2)                                        //
//----------------------------------------------------------------------//
void SecondaryForces(float fbb[2][Nb],float xb[2][Nb],int Nb,float hb,float Spr,int connect) {
  float dl1;
  float dl2;
  float dr1;
  float dr2;
  float ndl;
  float ndr;
  float Lrest;
  for (int ii=0; ii<Nb; ii++) {
    if (connect==1) {
      Lrest=2*hb;
      if (ii==0) {
        dl1=xb[0][Nb-2]-xb[0][ii];
        dl2=xb[1][Nb-2]-xb[1][ii];
        dr1=xb[0][ii+2]-xb[0][ii];
        dr2=xb[1][ii+2]-xb[1][ii];
      } else if (ii==1) {
        dl1=xb[0][Nb-1]-xb[0][ii];
        dl2=xb[1][Nb-1]-xb[1][ii];
        dr1=xb[0][ii+2]-xb[0][ii];
        dr2=xb[1][ii+2]-xb[1][ii];
      } else if (ii==Nb-2) {
        dl1=xb[0][ii-2]-xb[0][ii];
        dl2=xb[1][ii-2]-xb[1][ii];
        dr1=xb[0][0]-xb[0][ii];
        dr2=xb[1][0]-xb[1][ii];
      } else if (ii==Nb-1) {
        dl1=xb[0][ii-2]-xb[0][ii];
        dl2=xb[1][ii-2]-xb[1][ii];
        dr1=xb[0][1]-xb[0][ii];
        dr2=xb[1][1]-xb[1][ii];
      } else {
        dl1=xb[0][ii-2]-xb[0][ii];
        dl2=xb[1][ii-2]-xb[1][ii];
        dr1=xb[0][ii+2]-xb[0][ii];
        dr2=xb[1][ii+2]-xb[1][ii];
      }
      ndl=sqrt(pow(dl1,2)+pow(dl2,2));
      ndr=sqrt(pow(dr1,2)+pow(dr2,2));
      fbb[0][ii]=fbb[0][ii]+(Spr*(ndl-Lrest)*dl1/ndl)+(Spr*(ndr-Lrest)*dr1/ndr);
      fbb[1][ii]=fbb[1][ii]+(Spr*(ndl-Lrest)*dl2/ndl)+(Spr*(ndr-Lrest)*dr2/ndr);
    }
  }
  return;
} // SecondaryForces
//--------------------------------------------------------------------//

//----------------------------------------------------------------------//
// define central forces between boundary points and the cell nucleus   //
// only holf of all boundary points are connected to cell nucleus to    //
// allow the whole cell to grow                                         //
//----------------------------------------------------------------------//
void CenterForces(float fbb[2][Nb],float xb[2][Nb],int Nb,float cen,float len,float Spr,int connect){
  float dl1;
  float dl2;
  float ndl;
  float Lrest;

  for (int ii=0; ii<Nb/2; ii++){
    if (connect==2){
      Lrest=len;

      dl1=cen-xb[0][ii];
      dl2=cen-xb[1][ii];
      ndl=sqrt(pow(dl1,2)+pow(dl2,2));

      fbb[0][ii]=fbb[0][ii]+(Spr*(ndl-Lrest)*dl1/ndl);
      fbb[1][ii]=fbb[1][ii]+(Spr*(ndl-Lrest)*dl2/ndl);
    }
  }
  return;
} // function CenterForces
//--------------------------------------------------------------------//

//--------------------------------------------------------------------//
// define opposite forces between opposite boundary points            //
//--------------------------------------------------------------------//
void OppositeForces(float fbb[2][Nb],float xb[2][Nb],int Nb,float len,float Spr,int connect){
  float Lrest;
  float dl1;
  float dl2;
  float dr1;
  float dr2;
  float ndl;
  float ndr;
  int Nb2;

  for (int ii=0; ii<Nb; ii++){
    if (connect==3){
      Lrest=2*len;
      Nb2=floor(Nb/4);
      if (ii <= (Nb2+1)){
        dl1=xb[0][3*Nb2+2-ii]-xb[0][ii];
        dl2=xb[1][3*Nb2+2-ii]-xb[1][ii];
        dr1=-dl1;
        dr2=-dl2;
        ndl=sqrt(pow(dl1,2)+pow(dl2,2));
        ndr=sqrt(pow(dr1,2)+pow(dr2,2));
        fbb[0][ii]=fbb[0][ii]+(Spr*(ndl-Lrest)*dl1/ndl);
        fbb[1][ii]=fbb[1][ii]+(Spr*(ndl-Lrest)*dl2/ndl);
        fbb[0][3*Nb2+2-ii]=fbb[0][3*Nb2+2-ii]+(Spr*(ndr-Lrest)*dr1/ndr);
        fbb[1][3*Nb2+2-ii]=fbb[1][3*Nb2+2-ii]+(Spr*(ndr-Lrest)*dr2/ndr);
      }
    }
  }
  return;
}// function OppositeForces
//--------------------------------------------------------------------//


//-----------------------------------------------------------------//
// spreads the material values sb(Nb) (forces, sources) defined at //
// material points xb(1,Nb) to the fluid grid sg(Ng+1,Ng+1) in the //
// square domain [xmin,xmax]^2 with mesh width hg, material points //
// separation hb and a radius of the discrete delta function hdl.  //
//-----------------------------------------------------------------//
void BoundToGrid1(float sg[Ng+1][Ng+1],float xb[2][2],float sb[2][2],int Nb,int Ng,float hdl,float hg,float hb,int xmin,int xmax){
  int pas=-100; // passive value - do nothing
  float llx,rr,dx,lly,dy;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  int xbb1;
  int xbb2;
  int Nx,Ny;

  for (int n3=0;n3<Nb;n3++){
    // move points into the domain
    xbb1=IntoDom(xb[0][n3],xmin,xmax);
    xbb2=IntoDom(xb[1][n3],xmin,xmax);

    // determine indices of the nearest lower-down grid point
    Nx=1+floor((xbb1-xmin)/hg);
    Ny=1+floor((xbb2-xmin)/hg);

    // tests all 16 possible grid points
    for (int ii=-1; ii<2; ii++){
      for (int jj=-1; jj<2; jj++){
        // compute the interpolation Delta function
        llx=xmin+(Nx-1)*hg+ii*hg;
        rr=fabs(xbb1-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmin+(Ny-1)*hg+jj*hg;
        rr=fabs(xbb2-lly);
        dy=DeltaFun(rr,hdl);

        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);

        // update the values if poits are not pasive
        if (dx*dy > 0){
          sg[x1][y1]  += sb[0][n3]*dx*dy*hb;
          if (x2 != pas){
            sg[x2][y1]+= sb[0][n3]*dx*dy*hb;
          }
          if (y2 != pas){
            sg[x1][y2]+= sb[0][n3]*dx*dy*hb;
          }
          if ((x2 != pas) & (y2 != pas)){
            sg[x2][y2]+= sb[0][n3]*dx*dy*hb;
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
// square domain [xmin,xmax]^2 with mesh width hg, material points   //
// separation hb and a radius of the discrete delta function hdl.    //
//-------------------------------------------------------------------//
void BoundToGrid2(float sg[Ng+1][Ng+1][2],float xb[2][Nb],float sb[2][Nb],int Nb,int Ng,float hdl,float hg,float hb,int xmin,int xmax){

  // passive value - do nothing
  int pas=-100;
  int Nx,Ny,xbb1,xbb2;
  int x1=0;
  int x2=0;
  int y1=0;
  int y2=0;
  float llx,rr,dx,lly,dy;

  for (int n3=0; n3<Nb; n3++){
    // move points into the domain
    xbb1=IntoDom(xb[0][n3],xmin,xmax);
    xbb2=IntoDom(xb[1][n3],xmin,xmax);

    // determine indeces of the nearest lower-down grid point
    Nx=1+floor((xbb1-xmin)/hg);
    Ny=1+floor((xbb2-xmin)/hg);

    // tests all 16 possible grid points
    for (int ii=-1; ii<2;ii++){
      for (int jj=-1;jj<2;jj++){
        // compute the interpolation Delta function
        llx=xmin+(Nx-1)*hg+ii*hg;
        rr= fabs(xbb1-llx);
        dx= DeltaFun(rr,hdl);
        lly=xmin+(Ny-1)*hg+jj*hg;
        rr= fabs(xbb2-lly);
        dy= DeltaFun(rr,hdl);

        // determine indices of the grid points to update
        IndDel(x1,x2,llx,ii,Nx,Ng,xmin,xmax);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmin,xmax);

        // update the values if poits are not pasive
        if (dx*dy > 0){
          sg[x1][y1][1]=sg[x1][y1][1]+sb[0][n3]*dx*dy*hb;
          sg[x1][y1][2]=sg[x1][y1][2]+sb[1][n3]*dx*dy*hb;

          if (x2 != pas){
            sg[x2][y1][1]=sg[x2][y1][1]+sb[0][n3]*dx*dy*hb;
            sg[x2][y1][2]=sg[x2][y1][2]+sb[1][n3]*dx*dy*hb;
          }
          if (y2 != pas){
            sg[x1][y2][1]=sg[x1][y2][1]+sb[0][n3]*dx*dy*hb;
            sg[x1][y2][2]=sg[x1][y2][2]+sb[1][n3]*dx*dy*hb;
          }
          if ((x2 != pas) && (y2 != pas)){
            sg[x2][y2][1]=sg[x2][y2][1]+sb[0][n3]*dx*dy*hb;
            sg[x2][y2][2]=sg[x2][y2][2]+sb[1][n3]*dx*dy*hb;
          }
        }

      } // jj
    } // ii
  } // for n3
  return;
} // function BoundToGride
//---------------------------------------------------------------------//

//---------------------------------------------------------------------//
// Uses the Fast Fourier Method to solve the Navier-Stokes equations   //
// for updating grid velocity ug due to grid forces fg and grid source //
// distribution sg, rho and mu are fluid constants, dt is a time step, //
// hg is a mesh width.                                                 //
//---------------------------------------------------------------------//
void NavierStokes(float vg[Ng+1][Ng+1][2],float ug[Ng+1][Ng+1][2],float fg[Ng+1][Ng+1][2],float sg[Ng+1][Ng+1],int Ng,float rho,float mu,float dt,float hg){

  float pom,B1,B2,Aa,Bb,Bv;
  int in1,in2;
  float Eps=0.0000001;
  float sg_slice[Ng][Ng];
  float vg_slice1[Ng][Ng];
  float vg_slice2[Ng][Ng];
  float fvg1[Ng][Ng][2];
  float fvg2[Ng][Ng][2];
  float fvgg1[Ng][Ng];
  float fvgg2[Ng][Ng];
  float speqsg_slice[Ng][Ng/2];
  float speqvg_slice1[Ng][Ng/2];
  float speqvg_slice2[Ng][Ng/2];

  // stage n terms: force density fg, source distribution sg and current
  // velocity ug
  for (int n1=0; n1<Ng+1; n1++){
    for (int n2=0; n2<Ng+1; n2++){
      for (int ik=0; ik<2; ik++){
        // upwind scheme for the advection term
        if (ug[n1][n2][1] < 0){
          in1=PeriodInd(n1,Ng,1);
          pom=ug[in1][n2][ik]-ug[n1][n2][ik];
        }else{
          in1=PeriodInd(n1,Ng,-1);
          pom=ug[n1][n2][ik]-ug[in1][n2][ik];
        }
        vg[n1][n2][ik]=ug[n1][n2][1]*pom;

        if (ug[n1][n2][2] < 0){
          in2=PeriodInd(n2,Ng,1);
          pom=ug[n1][in2][ik]-ug[n1][n2][ik];
        }else{
          in2=PeriodInd(n2,Ng,-1);
          pom=ug[n1][n2][ik]-ug[n1][in2][ik];
        }
        vg[n1][n2][ik]=vg[n1][n2][ik]+ug[n1][n2][2]*pom;
        vg[n1][n2][ik]=-dt*vg[n1][n2][ik]/hg;

        // central difference for the grad of source term
        if (ik == 1){
          in1=PeriodInd(n1,Ng,1);
          in2=PeriodInd(n1,Ng,-1);
          pom=sg[in1][n2]-sg[in2][n2];
        }else if (ik == 2){
          in1=PeriodInd(n2,Ng,1);
          in2=PeriodInd(n2,Ng,-1);
          pom=sg[n1][in1]-sg[n1][in2];
        }
        vg[n1][n2][ik]=vg[n1][n2][ik]+dt*mu*pom/(6*hg*rho*rho);

        // current vlocity and force terms
        vg[n1][n2][ik]=vg[n1][n2][ik]+ug[n1][n2][ik]+dt*fg[n1][n2][ik]/rho;
      } // for ik
    } // for n2
  } // for n1

  // the Fast Fourier transforms of source distribution sg and stage n term vg
  for (int i=0;i<Ng;i++){
    for (int j=0;j<Ng;j++){
      sg_slice[i][j]  = sg[i][j];
      vg_slice1[i][j] = vg[i][j][1];
      vg_slice2[i][j] = vg[i][j][2];
    }
  }
  rlft3(sg_slice,speqsg_slice,1,Ng,Ng,1);
  rlft3(vg_slice1,speqvg_slice1,1,Ng,Ng,1);
  rlft3(vg_slice2,speqvg_slice2,1,Ng,Ng,1);

  // determines fug - the Fourier Transform of the velocity field at the stage n+1
  for (int n1=0;n1<Ng-1;n1++){
    for (int n2=0;n2<Ng-1;n2++){
      B1=sin(2*pi*n1/Ng);
      B2=sin(2*pi*n2/Ng);
      Bb=pow(B1,2)+pow(B2,2);
      Aa=1+4*mu*dt*(pow(sin(pi*n1/Ng),2)+pow(sin(pi*n2/Ng),2))/(rho*hg*hg);

      if (Bb < Eps){
        fvg1[n1+1][n2+1][0]=real(vg_slice1[n1+1][n2+1])/Aa;
        fvg1[n1+1][n2+1][1]=imag(vg_slice1[n1+1][n2+1])/Aa;
        fvg2[n1+1][n2+1][0]=real(vg_slice2[n1+1][n2+1])/Aa;
        fvg2[n1+1][n2+1][1]=imag(vg_slice2[n1+1][n2+1])/Aa;
      }else{
        Bv=B1*real(vg_slice1[n1+1][n2+1])+B2*real(vg_slice2[n1+1][n2+1]);
        fvg1[n1+1][n2+1][0]=(Bb*real(vg_slice1[n1+1][n2+1])-B1*Bv)/(Aa*Bb);
        fvg2[n1+1][n2+1][0]=(Bb*real(vg_slice2[n1+1][n2+1])-B2*Bv)/(Aa*Bb);

        Bv=B1*imag(vg_slice1[n1+1][n2+1])+B2*imag(vg_slice2[n1+1][n2+1]);
        fvg1[n1+1][n2+1][1]=(Bb*imag(vg_slice1[n1+1][n2+1])-B1*Bv)/(Aa*Bb);
        fvg2[n1+1][n2+1][1]=(Bb*imag(vg_slice2[n1+1][n2+1])-B2*Bv)/(Aa*Bb);

        fvg1[n1+1][n2+1][0]=fvg1[n1+1][n2+1][0]+hg*B1*imag(sg_slice[n1+1][n2+1])/(Bb*rho);
        fvg1[n1+1][n2+1][1]=fvg1[n1+1][n2+1][1]-hg*B1*real(sg_slice[n1+1][n2+1])/(Bb*rho);
        fvg2[n1+1][n2+1][0]=fvg2[n1+1][n2+1][0]+hg*B2*imag(sg_slice[n1+1][n2+1])/(Bb*rho);
        fvg2[n1+1][n2+1][1]=fvg2[n1+1][n2+1][1]-hg*B2*real(sg_slice[n1+1][n2+1])/(Bb*rho);
      }
    } // for n2
  } // for n1

  // the inverse Fast Fourier Method of fvg
  for (int i=0;i<Ng;i++){
    for (int j=0;j<Ng;j++){
      fvgg1[i][j]=fvg1[i][j][0]+i*fvg1[i][j][1];
      rlft3(fvgg1,fvgg1,1,Ng,Ng,-1);
    }
  }
  for (int i=0;i<Ng;i++){
    for (int j=0;j<Ng;j++){
      fvgg2[i][j]=fvg2[i][j][0]+i*fvg2[i][j][1];
      rlft3(fvgg2,fvgg2,1,Ng,Ng,-1);
    }
  }
  for (int ii=0; ii<Ng; ii++){
    for (int jj=0; jj<Ng; jj++){
      vg[ii][jj][0]=real(fvgg1[ii][jj]);
      vg[ii][jj][1]=real(fvgg2[ii][jj]);
    }
  }
  for (int ii=0; ii<Ng; ii++){
    vg[Ng][ii][0]=vg[0][ii][0];
    vg[ii][Ng][1]=vg[ii][0][1];
  }
} // function NavierStokes
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
// interpolates the grid values fg(Ng,Ng,2) (velocities) to the material //
// values fb(1,Nb) defined at the material points xb(1,Nb) in the square //
// domain [xmn,xmx]^2 with mesh width hg and a radius of the discrete    //
// Dirac delta hdl.                                                      //
//-----------------------------------------------------------------------//
void GridToBound(float fb[2][Nb],float xb[2][Nb],int Nb,float fg[Ng+1][Ng+1][2],int Ng,float hdl,float hg,int xmn,int xmx){

  float llx,rr,dx,lly,dy;
  int Nx,Ny,xbb1,xbb2;
  int x1 = 0;
  int x2 = 0;
  int y1 = 0;
  int y2 = 0;


  for (int n3=0; n3<Nb; n3++){
    // moves points into the domain
    xbb1=IntoDom(xb[0][n3],xmn,xmx);
    xbb2=IntoDom(xb[1][n3],xmn,xmx);

    //computes the indices of the nearest down-left grid point
    Nx=1+floor((xbb1-xmn)/hg);
    Ny=1+floor((xbb2-xmn)/hg);

    // test all 16 neighboring points
    for (int ii=-1;ii<2;ii++){
      for (int jj=-1;jj<2;jj++){
        // determine the value of the interpolation Delta-function
        llx=xmn+(Nx-1)*hg+ii*hg;
        rr=fabs(xbb1-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmn+(Ny-1)*hg+jj*hg;
        rr=fabs(xbb2-lly);
        dy=DeltaFun(rr,hdl);

        // determine the indices of grid points gaining positive impact
        IndDel(x1,x2,llx,ii,Nx,Ng,xmn,xmx);
        IndDel(y1,y2,lly,jj,Ny,Ng,xmn,xmx);

        // update the values if inside the impact domain
        if (dx*dy > 0){
          fb[0][n3]=fb[0][n3]+fg[x1][y1][1]*dx*dy*hg*hg;
          fb[1][n3]=fb[1][n3]+fg[x1][y1][2]*dx*dy*hg*hg;
        }
      } // for jj
    } // for ii
  }// for n3
  return;
} // function GridToBound
//-------------------------------------------------------------------//

void main() {

  int xmin=-1;
  int xmax=1;
  int cen=(xmax+xmin)/2;      // fluid domain (square)
  float hg=(xmax-xmin)/Ng;    // fluid mesh width
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
  int NumLoop=40;             // number of steps
  int mod_num=5;              // frequency
  int Nbs=2;
  float sb[2][2];
  float xb[2][Nb] = {0};      // Position of all boundary elements
  float ub[2][Nb] = {0};      // Velocity of all boundary elements
  float hb;                   // Spacing between boundary elements?
  float fb[2][Nb] = {0};      // Boundary forces
  float fadj[2][Nb]={0};      // Adjacent forces array
  float fsec[2][Nb]={0};      // Secondary forces array
  float fcen[2][Nb]={0};      // Centre forces array
  float fopp[2][Nb]={0};      // Opposite forces array
  float xg[Ng][Ng][2];      // Fluid grid array
  float sg[Ng+1][Ng+1]={0};// Fluid grid
  float fg[Ng+1][Ng+1][2]={0};// Grid forces
  float vg[Ng+1][Ng+1][2]={0};// Grid velocities
  float ug[Ng+1][Ng+1][2]={0};// Previous grid velocities
  float sbb[2][2]  ={0};      //

  //-- define fluid grid --//
  for (int ii=0;ii<Ng;ii++) {
    for (int jj=0;jj<Ng;jj++) {
      xg[ii][jj][1]=xmin+(ii-1)*hg;
      xg[ii][jj][2]=xmin+(jj-1)*hg;
    }
  }

  //-- cell shape --//
  for (int ii=0;ii<Nb;ii++) {
    xb[1][ii] = DefineCellShape(Nb,len,1,ii);
    xb[2][ii] = DefineCellShape(Nb,len,2,ii);
  }
  hb=sqrt(pow((xb[1][1]-xb[1][2]),2)+pow((xb[2][1]-xb[2][2]),2));

  //-- distributed sources and sinks --//
  sb[0][0]=cen;
  sb[1][0]=cen;
  sb[0][1]=xmin;
  sb[1][1]=xmin;
  sbb[0][0]= Src;     // a source at the center
  sbb[0][1]=-Src;     // a sink in the corner

  for (int loop_num=0; loop_num<NumLoop; loop_num++) {
    //-- new position of boundary points --//
    for (int ii=0; ii<Nb; ii++){
      xb[0][ii] = xb[0][ii]+dt*ub[0][ii];
      xb[1][ii] = xb[1][ii]+dt*ub[1][ii];
    }
    //-- boundary forces --//
    AdjacentForces(fadj,xb,Nb,hb,Spr);            // adjacent
    SecondaryForces(fsec,xb,Nb,hb,Spr,connect);   // secondary
    CenterForces(fcen,xb,Nb,cen,len,Spr,connect); // center
    OppositeForces(fopp,xb,Nb,len,Spr,connect);   // opposite
    for (int ii=0; ii<Nb; ii++){
      fb[0][ii]=fadj[0][ii]+fsec[0][ii]+fcen[0][ii]+fopp[0][ii];    // add all forces
      fb[1][ii]=fadj[1][ii]+fsec[1][ii]+fcen[1][ii]+fopp[1][ii];    // add all forces
    }
    //-- grid sources --//
    BoundToGrid1(sg,sb,sbb,Nbs,Ng,hg,hg,0.5*hg,xmin,xmax);
    //-- grid forces --//

    BoundToGrid2(fg,xb,fb,Nb,Ng,hg,hg,0.5*hg,xmin,xmax);
    //-- compute grid velocity from NavierStokes --//
    NavierStokes(vg,ug,fg,sg,Ng,rho,mu,dt,hg);
    for (int ii=0;ii<Ng;ii++) {
      for (int jj=0;jj<Ng;jj++) {
        ug[ii][jj][1]=vg[ii][jj][1];
        ug[ii][jj][2]=vg[ii][jj][2];
      }
    }
    //-- boundary velocities --//
    GridToBound(ub,xb,Nb,vg,Ng,hg,hg,xmin,xmax);
  }   // for loop_num
}
