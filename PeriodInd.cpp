#include "routines.hpp"

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
