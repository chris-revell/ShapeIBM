#include <math.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

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
//Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3] (where nn1 = 1 for the case of a logically two-dimensional array), this routine returns (for isign=1) the complex fast Fourier transform as two complex arrays: On output, data contains the zero and positive frequency values of the third frequency component, while speq[1..nn1][1..2*nn2] contains the Nyquist critical frequency values of the third frequency component. First (and second) frequency components are stored for zero, positive,
//and negative frequencies, in standard wrap- around order. See text for description of how complex values are arranged. For isign=-1, the inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative factor) is performed, with output data (viewed as a real array) deriving from input data (viewed as complex) and speq. For inverse transforms on data not generated first by a forward transform, make sure the complex input data array satisfies property (12.5.2). The dimensions nn1, nn2, nn3 must always be integer powers of 2.
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

int main(){
  int i,j,k;
  int ndim1=1;
  int ndim2=64;
  int ndim3=64;

  // memory allocation
  float*** data=new float**[ndim1];
  for(i=0;i<ndim1;i++){
  	data[i]=new float*[ndim2];
  	for(j=0;j<ndim2;j++){
  		data[i][j]=new float[ndim3];
  	}
  }

  // access to 3D array
  for(i=0;i<ndim1;i++){
  	for(j=0;j<ndim2;j++){
  		for(k=0;k<ndim3;k++){
  			data[i][j][k]=rand();
  		}
  	}
  }

  float** speq=new float*[ndim1];
  for(i=0;i<ndim1;i++){
  	data[i]=new float*[ndim2];  	
  }

  // access to 3D array
  for(i=0;i<ndim1;i++){
  	for(j=0;j<ndim2;j++){
  		for(k=0;k<ndim3;k++){
  			data[i][j][k]=0;
  		}
  	}
  }

  rlft3(data,spec,1,4,4,1);


  // Freeing the memory
  for(i=0;i<ndim1;i++){
  	for(j=0;j<ndim2;j++){
  		delete[] data[i][j];
  	}
  	delete[] data[i];
  }
  delete[] data;

}
