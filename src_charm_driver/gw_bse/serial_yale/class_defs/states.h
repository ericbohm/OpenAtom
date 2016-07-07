// class STATES
// STATES has index for spin and k point.
// STATES class is defined as 2 dimensional array in main program

#ifndef STATES_H
#define STATES_H


#include "../include/ckcomplex.h"
#include "gspace.h"

class STATES{

  public:

    int ikpt;         // k index
    int ispin;        // spin index
    
    int nstate;       // number of states
    int nocc;         // number of occupied states
    int nunocc;       // number of unoccupied states

    int ndata;        // number of data (# of g vector or # of fft grid)
    
    bool shifted;     // if wavefunction shifted or not
    
    double *eig;      // eigenvalues (in Hartree)
    double *occ;      // occupation ( usually 1 (occupied) or 0 (unoccupied) )

    complex **coeff;  // wavefunction coefficient defined as 2 dimensional array
                      // coeff[istate][i] : first index is for state index
                      //                    second index is for g or fftgrid index
    
    GSPACE gvec;

    
    // constructor
    STATES(int _ispin, int _ikpt, int _nstate, bool _shifted){
	ispin = _ispin;
	ikpt = _ikpt;
	nstate = _nstate;
	shifted = _shifted;
	eig = new double [nstate];
	occ = new double [nstate];
    }
       
};


#endif
