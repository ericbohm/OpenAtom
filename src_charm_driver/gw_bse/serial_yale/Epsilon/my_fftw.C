#include <cstdlib>
#include <iostream>
#include "my_fftw.h"
#include "util.h"

// if we don't declare this, the code complains and don't compile
fftw_complex *in_pointer;
fftw_complex *out_pointer;
fftw_plan my_plan;

// execute fftw3
void do_fftw(){
    fftw_execute(my_plan);
}

// not so sure what the purpose of these two routines
fftw_complex *fftw_in_pointer(){
    return in_pointer;
}
fftw_complex *fftw_out_pointer(){
    return out_pointer;
}


// destroy fftw related variables. We only do once when fftw is all done
void destroy_fftw_stuff(){
    mymessage("fftw_plan and pointers are being destroyed");
    fftw_destroy_plan(my_plan);
    fftw_free(in_pointer);
    fftw_free(out_pointer);
}


// setup 3D fftw plan
void setup_fftw_3d(int nfft[3], int direction){

    static bool firsttime = true;
    static int nfft0_old = 0, nfft1_old = 0, nfft2_old = 0, direction_old = 0;

    // check for some dumb input values
    if (nfft[0]<=0 || nfft[1] <=0 || nfft[2] <=0){
	Die("setup_fftw_3d routine received illigal value for nfft. nfft should be positive number.");
    }
    if (!(direction==-1 || direction==1)){
	Die("setup_fftw_3d routine received illigal value for direction. FFTW direction must either 1 or -1");
    }

    // if first time, we need to set up
    if(firsttime){
	const int ndata = nfft[0]*nfft[1]*nfft[2];
	in_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
	out_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
	my_plan = fftw_plan_dft_3d(nfft[0], nfft[1], nfft[2], in_pointer, out_pointer, direction, FFTW_ESTIMATE);
    }


    // if not the first time and any parameters (size, direction) mismatch
    // we need to destroy old plans to set up new ones
    if (!firsttime && ((direction_old != direction) ||
		       (nfft[0] != nfft0_old) ||
		       (nfft[1] != nfft1_old) ||
		       (nfft[2] != nfft2_old) ) ){
	destroy_fftw_stuff();

	const int ndata = nfft[0]*nfft[1]*nfft[2];
	in_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
	out_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
	my_plan = fftw_plan_dft_3d(nfft[0], nfft[1], nfft[2], in_pointer, out_pointer, direction, FFTW_ESTIMATE);
    }

    // if not first time but perfect match of parameters, then we don't need to set up
    // since we did it already! Just do nothing.
    if(!firsttime && direction_old==direction && nfft0_old == nfft[0] &&
	  nfft1_old == nfft[1] && nfft2_old == nfft[2]){
    }

    // now the old value changes to new value
    firsttime = false;
    nfft0_old = nfft[0];
    nfft1_old = nfft[1];
    nfft2_old = nfft[2];
    direction_old = direction;

}
