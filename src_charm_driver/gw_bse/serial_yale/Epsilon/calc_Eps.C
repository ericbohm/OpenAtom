/* This routine calculates Epsilon matrix
   Input: Polarizability matrix
   nrow/ncol in Epsmat are smaller (about half) than nrow/ncol in P since we do spherical cutoff for Epsilon
   bool *accept tells us which G vectors we are going to use for Epsilon matrix
*/


#include <iostream>
#include <cstdlib>
#include "include/ckcomplex.h"
#include "class_defs/matrix.h"
#include "coulomb.h"
#include "util.h"


void calc_Epsmat(double* vcoulb, int iq, SYSINFO sys, GSPACE *geps, CMATRIX *Epsmat, CMATRIX *P, int ndata, bool* accept){

    vcoulb = new double [geps->ng];
    
    // calculate Vcoulb
    calc_vcoulb(iq, vcoulb, geps, sys);
    
    int inew=0;
    int jnew=0;
    
    for (int i=0; i<ndata; i++) {
        if (accept[i]) {
            jnew = 0;
            for (int j=0; j<ndata; j++) {
                if (accept[j]) {
                    Epsmat->get(inew,jnew) = P->get(i,j);
		    
                    Epsmat->get(inew,jnew) *= -1*sqrt( vcoulb[inew] )*sqrt( vcoulb[jnew] );
                    if ( inew == jnew ){
                        Epsmat->get(inew,jnew) += double(1);
                    }
                    jnew += 1;
                }
            }
        inew += 1;
        }
    }

    // print out errors, probably die because it would be a critical error
    if (inew!=geps->ng) {
        Die(" Number of G doesn't match with the number of rows in Epsmat!!");
    }
    else if (jnew!=geps->ng){
        Die(" Number of G doesn't match with the number of columns in Epsmat!!");
    }
    
    delete[] vcoulb;
    
}
