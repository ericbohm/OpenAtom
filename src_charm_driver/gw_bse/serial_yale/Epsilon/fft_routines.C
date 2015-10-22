#include "fft_routines.h"
/* fft subroutines */

// g index to fft index
// This function takes g list from the wavefunction (in momentum space) and
// get fft index
// We do this because ng is ALWAYS less than nfft[0]*nfft[1]*nfft[2]
void gidx_to_fftidx(int ng, int **g, int nfft[3], int **fftidx){
    
    int this_gidx[3];
    int keep_gidx[3];
    int upper, lower, sumi;

    // loop over g index
    for(int ig=0; ig<ng; ig++){
	
        // g index at ig
        this_gidx[0] = g[ig][0];
        this_gidx[1] = g[ig][1];
        this_gidx[2] = g[ig][2];
        
        // speculate if this_gidx[3] is inside of the fftbox
	// keep_gidx[i] has value -1 if this_gidx[i] is outside of the fftbox
	// keep_gidx[i] has value 0 if this_gidx[i] is inside of the fftbox
	// here, we initialize keep_gidx with -1
        for (int i=0;i<3;i++){ keep_gidx[i]=-1; }

        sumi = 0;
        for (int j=0; j<3; j++) {
	    // upper bound of the fftbox
            upper = nfft[j]/2;
	    // lower bound of the fftbox
            lower = -1*nfft[j]/2;
            if (this_gidx[j]<upper && this_gidx[j] >= lower) {
                keep_gidx[j] = 0;
            }
            sumi += keep_gidx[j];
        }

	// we put g index to fftbox index if keep_gidx[i] is 0 for all i
        if (sumi == 0) {
            for (int j=0; j<3; j++) {
                fftidx[ig][j] = this_gidx[j] + 1;
		// we don't want our fftbox index becomes negative
                // make all fftidx positive
                if (fftidx[ig][j]<=0) {
                    fftidx[ig][j] += nfft[j];
                }
            }
        }
        else{
            // we need this to reject g if this_gidx doesn't fit into fftbox
            for (int i=0;i<3;i++){ fftidx[ig][i] = -1; }
        }
        
    }// end ig loop
    
}




// this routine takes g space wavefunction (array) and assign values into 1D array 
void put_into_fftbox(int ng, complex *rawdata, int **fftidx, int nfft[3], fftw_complex *in_pointer){

    // initialize in
    int ndata=nfft[0]*nfft[1]*nfft[2];
    for (int i=0; i<ndata; i++){
	in_pointer[i][0]=0;
	in_pointer[i][1]=0;
    }
    
    // temporary index
    int idxtmp;
    for (int ig=0; ig<ng; ig++){
	idxtmp = (fftidx[ig][0]-1)*nfft[1]*nfft[2] +
                    (fftidx[ig][1]-1)*nfft[2] +
                     fftidx[ig][2];
        // assign wavefunction value to in
	// we subtract 1 since our arrays use C++ counting
        in_pointer[idxtmp-1][0] = rawdata[ig].re;
	in_pointer[idxtmp-1][1] = rawdata[ig].im;
    }
}



// after fftw_execute, we want to transfer data from "out" to "array"
// you can scale the values because fftw doesn't do any scale for you
void fftbox_to_array(int ndata, fftw_complex *out_pointer, complex *array, double scale){
    for(int i=0; i<ndata; i++){
	array[i].re = out_pointer[i][0]*scale;
	array[i].im = out_pointer[i][1]*scale;
    }
    
}


// fft index is all positive integers, so this routine changes fftindex to G index
void fftidx_to_gidx(int *gx, int *gy, int *gz, int nfft[3]){
    
    int ijk=0;
    
    for (int i=0; i<nfft[0]; i++) {
        for (int j=0; j<nfft[1]; j++) {
            for (int k=0; k<nfft[2]; k++) {
                gx[ijk] = i;
                gy[ijk] = j;
                gz[ijk] = k;
                ijk += 1;
            }
        }
    }
   
   int ndata = nfft[0]*nfft[1]*nfft[2];
    
    for (int ijk=0; ijk<ndata; ijk++) {
        if (gx[ijk] >= nfft[0]/2) {
            gx[ijk] -= nfft[0];
        }
        if (gy[ijk] >= nfft[1]/2) {
            gy[ijk] -= nfft[1];
        }
        if (gz[ijk] >= nfft[2]/2) {
            gz[ijk] -= nfft[2];
        }
    }
    
    
}
