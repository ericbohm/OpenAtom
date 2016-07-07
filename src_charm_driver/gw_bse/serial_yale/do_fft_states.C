// do_fft_states performs fft for all k points and states (bands)
// For each spin channel, do_fft_states has to be called 

#include <cstdlib>
#include <iostream>
#include "class_defs/sysinfo.h"
#include "class_defs/states.h"
#include "my_fftw.h"
#include "fft_routines.h"
#include "util.h"


void do_fft_states(STATES **psi, STATES **psiR, SYSINFO sys, int nfft[3]){

    const int nkpt = sys.nkpt;
    const int nstate = psi[0]->nstate;

    // number of data (=number of fft grids)
    const int ndata = nfft[0]*nfft[1]*nfft[2];

    // number of g vectors
    int ng;
    // g vectors (double pointer)
    int **g;
    
    // scale
    double scale = 1/sqrt(ndata);
   
    setup_fftw_3d(nfft,1);
    
    for (int ik=0; ik < nkpt; ik++){
	// get g list     
        // TODO : this can go outside of this k loop since we now have the same number of g in each k point
	// 1. number of g vectors
	ng = psi[ik]->gvec.ng;

	// 2. m'alloc g vectors, and get data from psi[is][ik][0]
	// (remember! g vectors are saved only at the first state)
	g = new int *[ng];
	for (int ig=0; ig < ng; ig++){
	    g[ig] = new int [3];
	    g[ig][0] = psi[ik]->gvec.ig[ig];
	    g[ig][1] = psi[ik]->gvec.jg[ig];
	    g[ig][2] = psi[ik]->gvec.kg[ig];
	}

	
	// m'alloc fftindex 
	
	int **fftidx;
	fftidx = new int *[ng];
	for (int ig=0; ig < ng; ig++){
	    fftidx[ig] = new int [3];
	}
	

	
	// 3. transfer g index to fftindex
	gidx_to_fftidx(ng, g, nfft, fftidx); 
	// m'alloc psiR wavefunction coefficient
	psiR[ik]->coeff = new complex *[nstate];

	for(int ib=0; ib < nstate; ib++){

	    put_into_fftbox(ng, psi[ik]->coeff[ib], fftidx, nfft, in_pointer);

	    // call fftw
	    do_fftw();

	    // m'alloc psiR wavefunction
	    psiR[ik]->coeff[ib] = new complex [ndata];

	    fftbox_to_array(ndata, out_pointer, psiR[ik]->coeff[ib], scale);

	}
    }


    //destroy_fftw_stuff();

    mymessage("FFT has been performed");

    // let's copy eigenvalues and occupancies here
    for(int ik=0; ik<nkpt; ik++){
	for(int ib=0; ib<nstate; ib++){
	    psiR[ik]->eig[ib] = psi[ik]->eig[ib];
	    psiR[ik]->occ[ib] = psi[ik]->occ[ib];
	}
        psiR[ik]->ndata = ndata;
    }

}
