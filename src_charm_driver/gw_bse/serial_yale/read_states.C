/*
      read_states.C
        1) read eignevalues and occupancies from file Eigenvalue
        1) read wavefunctions from state*.out files
        2) read gvectors from state1.out file. G vectors are saved only once when state index is 0 at each k point
*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include "class_defs/states.h"
#include "util.h"

using namespace std;

void read_states(STATES *psi){

    const int nstate = psi->nstate;

    // memory allocation for wavefunction coefficient
    psi->coeff = new complex *[nstate];
    
    // string to open files 
    stringstream ss; // spin index
    stringstream sk; // kpoint index
    stringstream sb; // state index
    string fname;    // file name
    ifstream fp;     // file pointer

    ss << psi->ispin;
    sk << psi->ikpt;

    // number of g vectors
    int ng;
    // number of dense FFT grid points for density
    int nfftDen[3];
    
    // loop over states
    for (int istate=0; istate < nstate; istate++){

	// state
	sb << istate + 1;

	// if the state is not shifted
	if (!psi->shifted)
	    fname =  "./STATES_IN/Spin." + ss.str() + "_Kpt." + sk.str() + "_Bead.0_Temper.0/state" + sb.str() + ".out";
	// if the state is not shifted
	if (psi->shifted)
	    fname = "./STATES_IN/Spin." + ss.str() + "_Kpt.0" + sk.str() + "_Bead.0_Temper.0/state" + sb.str() + ".out";

	fp.open( fname.c_str() );

	// read ng and dense fft size
	fp >> ng >> nfftDen[0] >> nfftDen[1] >> nfftDen[2];
	if (istate == 0) {
          psi->ndata = ng;
        }
	// m'alloc for coefficient
	psi->coeff[istate] = new complex [ng];

	// read wavefunction and gvectors

	// g vectors are read when istate=0
	if (istate==0){
	    psi->gvec.ng = ng;
	    psi->gvec.ig = new int [ng];
	    psi->gvec.jg = new int [ng];
	    psi->gvec.kg = new int [ng];
	    for (int i=0; i<ng; i++){
		fp >> psi->coeff[istate][i].re >> psi->coeff[istate][i].im >> psi->gvec.ig[i] >> psi->gvec.jg[i] >> psi->gvec.kg[i];
	    }
	}
	else{
	    int dummy;
	    for (int i=0; i<ng; i++){
		fp >> psi->coeff[istate][i].re >> psi->coeff[istate][i].im >> dummy >> dummy >> dummy;
	    }
	}
	// close state file
	fp.close();

	// clean the stringstream
	sb.str("");
    }// end istate loop


    // read eigenvalues and occupancies

    if(!psi->shifted)
	fname = "./STATES_IN/Spin." + ss.str() + "_Kpt." + sk.str() + "_Bead.0_Temper.0/eigenvalues.in";
    if(psi->shifted)
        fname = "./STATES_IN/Spin." + ss.str() + "_Kpt.0" + sk.str() + "_Bead.0_Temper.0/eigenvalues.in";

    fp.open( fname.c_str() );
    
    if(!fp)
	Die("Failed to open Eigenvalue file");
    
    
    if(fp){
        for (int i=0; i < nstate; i++) {
            fp >> psi->eig[i];
        }
    }
    fp.close();
    
}

